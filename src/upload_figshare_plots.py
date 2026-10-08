#!/usr/bin/env python3
"""Resumable, rate-limited Figshare draft uploader. Never publishes records."""
import argparse
import csv
import fcntl
import hashlib
import html
import json
import sqlite3
import threading
import time
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
from pathlib import Path
from urllib.parse import urlparse

import requests

BASE = "https://api.figshare.com/v2/"


class VerificationPending(RuntimeError):
    """The transfer completed but Figshare is still processing it."""


class Client:
    def __init__(self, token):
        self.token = token
        self.local = threading.local()
        self.rate_lock = threading.Lock()
        self.last = 0

    def request(self, method, endpoint, **kwargs):
        url = endpoint if endpoint.startswith("https://") else BASE + endpoint
        host = urlparse(url).hostname
        if urlparse(url).scheme != "https" or not (host == "figshare.com" or host.endswith(".figshare.com")):
            raise RuntimeError("Unexpected API/upload host; refusing request")
        headers = {"Authorization": "token " + self.token} if host == "api.figshare.com" else {}
        if not hasattr(self.local, "session"):
            self.local.session = requests.Session()
        # Only safe/idempotent requests retry automatically. On uncertain POST
        # results, stop; restarting reconciles titles and files with the server.
        for attempt in range(6):
            with self.rate_lock:
                time.sleep(max(0, 1.05 - (time.monotonic() - self.last)))
                self.last = time.monotonic()
            try:
                response = self.local.session.request(method, url, headers=headers,
                                                timeout=(30, 180), allow_redirects=False, **kwargs)
            except requests.RequestException:
                if method not in ("GET", "PUT") or attempt == 5:
                    raise RuntimeError(f"Network failure during {method}; safe to resume") from None
                time.sleep(min(60, 2 ** (attempt + 1)))
                continue
            if response.status_code == 429 or (response.status_code >= 500 and method in ("GET", "PUT")):
                if attempt < 5:
                    time.sleep(min(60, max(2 ** (attempt + 1), int(response.headers.get("Retry-After", "0")))))
                    continue
            if not 200 <= response.status_code < 300:
                # Do not print response bodies or signed upload URLs.
                route = urlparse(url).path if host == "api.figshare.com" else "upload-service"
                raise RuntimeError(f"Figshare HTTP {response.status_code} during {method} {route}; safe to resume")
            if not response.content:
                return None
            try:
                return response.json()
            except ValueError:
                if (method == "PUT" and host != "api.figshare.com") or (
                    method == "POST" and "/files/" in urlparse(url).path
                ):
                    return None
                raise RuntimeError(f"Unexpected non-JSON response: {method}, host {host}, HTTP {response.status_code}") from None
        raise RuntimeError("Request retries exhausted")


def digest(path):
    h = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def complete(info, md5):
    return info.get("computed_md5") == md5 and (
        info.get("status") == "available" or info.get("is_incomplete") is False
    )


def transfer_file(client, staging, endpoint, row, info):
    """Transfer only; Figshare's asynchronous checks never block this worker."""
    path = staging / row["staged_path"]
    if path.stat().st_size != int(row["bytes"]):
        raise RuntimeError("Source size changed since staging")
    md5 = digest(path)
    db = sqlite3.connect(staging / "upload.sqlite", timeout=60)
    try:
        saved = db.execute("SELECT id, md5 FROM files WHERE path=?", (row["staged_path"],)).fetchone()
        if info and info.get("supplied_md5") != md5:
            raise RuntimeError("Existing file checksum differs; refusing overwrite")
        transferred = info and saved and saved == (info["id"], md5)
        processing = info and info.get("status") in ("ic_checking", "available")
        new_transfer = False
        if not (info and (complete(info, md5) or transferred or processing)):
            if info:
                info = client.request("GET", f"{endpoint}/files/{info['id']}")
            else:
                result = client.request("POST", endpoint + "/files", json={
                    "name": row["filename"], "size": int(row["bytes"]), "md5": md5})
                info = client.request("GET", result["location"])
            if not complete(info, md5) and info.get("status") not in ("ic_checking", "available"):
                upload_url = info["upload_url"]
                parts = client.request("GET", upload_url)
                with path.open("rb") as handle:
                    for part in parts["parts"]:
                        if part["status"] == "COMPLETE":
                            continue
                        handle.seek(part["startOffset"])
                        payload = handle.read(part["endOffset"] - part["startOffset"] + 1)
                        client.request("PUT", upload_url + "/" + str(part["partNo"]), data=payload)
                client.request("POST", f"{endpoint}/files/{info['id']}")
                new_transfer = True
        if info.get("status") == "available" and info.get("computed_md5") and info["computed_md5"] != md5:
            raise RuntimeError(f"Checksum mismatch for file {info['id']}")
        db.execute("INSERT OR REPLACE INTO files VALUES (?, ?, ?, ?)",
                   (row["staged_path"], info["id"], md5, int(complete(info, md5))))
        db.commit()
        return new_transfer
    finally:
        db.close()


def verify_record(client, db, endpoint, expected, require_all=True):
    """One metadata request verifies all available files in a record."""
    article = client.request("GET", endpoint)
    files = {f["name"]: f for f in article["files"]}
    expected_names = {r["filename"] for r in expected}
    if set(files) - expected_names:
        raise RuntimeError("Unexpected remote files during verification")
    pending = []
    for row in expected:
        saved = db.execute("SELECT id, md5 FROM files WHERE path=?", (row["staged_path"],)).fetchone()
        if not saved:
            if require_all:
                raise RuntimeError("Missing transfer checkpoint")
            continue
        info = files.get(row["filename"])
        if not info or info["id"] != saved[0]:
            raise RuntimeError("Remote file missing or replaced during verification")
        if info.get("status") == "available" and info.get("computed_md5") and info["computed_md5"] != saved[1]:
            raise RuntimeError(f"Checksum mismatch for file {info['id']}")
        verified = complete(info, saved[1])
        db.execute("UPDATE files SET verified=? WHERE path=?", (int(verified), row["staged_path"]))
        if not verified:
            pending.append({"id": info["id"], "name": row["filename"], "status": info.get("status", "unknown")})
    db.commit()
    return pending


def parallel_transfers(transfer, rows, workers):
    """Keep at most workers transfers in flight; stop scheduling on errors."""
    iterator = iter(rows)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        active = set()
        for _ in range(workers):
            row = next(iterator, None)
            if row is not None:
                active.add(pool.submit(transfer, row))
        while active:
            done, active = wait(active, return_when=FIRST_COMPLETED)
            results = [future.result() for future in done]
            for result in results:
                yield result
                row = next(iterator, None)
                if row is not None:
                    active.add(pool.submit(transfer, row))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--token-file", type=Path, required=True)
    parser.add_argument("--staging", type=Path, default=Path("int/figshare_pooled_upload"))
    parser.add_argument("--limit", type=int, help="Maximum new uploads, for a pilot")
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    if args.workers < 1 or args.workers > 8:
        parser.error("--workers must be between 1 and 8")
    staging = args.staging.resolve()
    if (staging / "SUPERSEDED").exists():
        raise RuntimeError("This upload manifest was superseded; use pooled-only staging")
    lock = (staging / "upload.lock").open("w")
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    if args.token_file.stat().st_mode & 0o077:
        raise RuntimeError("Token file must have owner-only permissions")
    client = Client(args.token_file.read_text().strip())
    records = list(csv.DictReader((staging / "records.csv").open()))
    rows = list(csv.DictReader((staging / "manifest.csv").open()))
    db = sqlite3.connect(staging / "upload.sqlite")
    db.execute("CREATE TABLE IF NOT EXISTS records (name TEXT PRIMARY KEY, id INTEGER NOT NULL)")
    db.execute("CREATE TABLE IF NOT EXISTS files (path TEXT PRIMARY KEY, id INTEGER, md5 TEXT, verified INTEGER)")
    db.commit()
    account = client.request("GET", "account")
    total_bytes = sum(int(row["bytes"]) for row in rows)
    if total_bytes > account["quota"]:
        raise RuntimeError("Dataset exceeds account quota")
    print(f"Account verified; dataset {total_bytes:,} bytes; quota {account['quota']:,} bytes", flush=True)
    existing = []
    for page in range(1, 7):
        batch = client.request("GET", f"account/articles?page_size=100&page={page}")
        existing.extend(batch)
        if len(batch) < 100:
            break
    by_title = {}
    for article in existing:
        by_title.setdefault(article["title"], []).append(article)
    links = []
    for record in records:
        title = "Synapse scatter plots — " + record["title"]
        saved = db.execute("SELECT id FROM records WHERE name=?", (record["record"],)).fetchone()
        if saved:
            article_id = saved[0]
        elif title in by_title:
            if len(by_title[title]) != 1:
                raise RuntimeError("Duplicate matching draft titles; manual reconciliation required")
            article_id = by_title[title][0]["id"]
        else:
            description = (
                f"<p>Individual synapse scatter and depth-density plots: {html.escape(record['title'])}.</p>"
                "<p>Filenames identify gene, anatomical hemisphere, presynaptic or postsynaptic "
                "compartment, and Notch condition. Only pooled synapse density is included, using "
                "pooled synapse depths. "
                "Legend PDFs are included. Original PDF contents are unchanged; files are not merged.</p>"
            )
            result = client.request("POST", "account/articles", json={"title": title, "description": description})
            article_id = int(result["location"].rstrip("/").split("/")[-1])
        db.execute("INSERT OR REPLACE INTO records VALUES (?, ?)", (record["record"], article_id))
        db.commit()
        links.append(dict(record=record["record"], article_id=article_id,
                          edit_url=f"https://figshare.com/account/articles/{article_id}/edit", files=record["files"]))
        # Write incrementally so links are available before uploads finish.
        with (staging / "draft_links.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(links[0]))
            writer.writeheader(); writer.writerows(links)
        print(f"Draft {article_id}: {record['title']}", flush=True)
    if args.prepare_only:
        return
    count = 0
    for link in links:
        article_id = link["article_id"]
        endpoint = f"account/articles/{article_id}"
        article = client.request("GET", endpoint)
        if article.get("published_date") or article.get("is_public"):
            raise RuntimeError("Refusing to modify a published record")
        remote = {f["name"]: f for f in article["files"]}
        expected = [row for row in rows if row["record"] == link["record"]]
        expected_names = {row["filename"] for row in expected}
        if set(remote) - expected_names:
            raise RuntimeError("Unexpected files found in draft; refusing to alter it")
        work = expected
        if args.limit:
            saved_paths = {r[0] for r in db.execute("SELECT path FROM files")}
            work = [r for r in expected if r["staged_path"] not in saved_paths][:max(0, args.limit-count)]
        # No per-file checksum polling: overlap transfers and server processing.
        def transfer(row):
            return transfer_file(client, staging, endpoint, row, remote.get(row["filename"]))
        for new_transfer in parallel_transfers(transfer, work, args.workers):
            count += int(new_transfer)
            if new_transfer and (count == 1 or count % 25 == 0):
                transferred = db.execute("SELECT COUNT(*) FROM files").fetchone()[0]
                print(f"Transferred {transferred}/{len(rows)} PDFs; {count} new this run", flush=True)
        pending = verify_record(client, db, endpoint, expected, require_all=not args.limit)
        print(f"Record {article_id}: transfer pass finished; {len(pending)} awaiting Figshare processing", flush=True)
        if args.limit and count >= args.limit:
            print("Pilot transfer pass complete; checks are batched, pending files do not block uploads", flush=True)
            return
    # Recheck only affected records after ALL transfers have been attempted.
    deadline = time.monotonic() + 1800
    remaining_links = links
    while remaining_links:
        pending_files = []
        retry_links = []
        for link in remaining_links:
            expected = [row for row in rows if row["record"] == link["record"]]
            pending = verify_record(client, db, f"account/articles/{link['article_id']}", expected)
            if pending:
                retry_links.append(link)
                pending_files.extend(pending)
        (staging / "pending_verification.json").write_text(json.dumps(pending_files, indent=2) + "\n")
        if not retry_links:
            break
        print(f"All transfer passes finished; {len(pending_files)} PDFs still processing", flush=True)
        if time.monotonic() >= deadline:
            raise VerificationPending(f"All transfers finished; {len(pending_files)} PDFs await Figshare checks")
        time.sleep(min(60, max(0, deadline-time.monotonic())))
        remaining_links = retry_links
    print(f"SUCCESS: {len(rows)} PDFs checksum verified in {len(records)} unpublished draft records", flush=True)
    (staging / "UPLOAD_COMPLETE.json").write_text(json.dumps({"files": len(rows), "records": len(records), "published": False}) + "\n")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        # Never expose request URLs, response bodies, credentials, or traceback locals.
        if isinstance(exc, RuntimeError):
            print(f"STOPPED: {exc}", flush=True)
        else:
            print(f"STOPPED: {type(exc).__name__}; inspect local inputs and resume", flush=True)
        raise SystemExit(75 if isinstance(exc, VerificationPending) else 1)
