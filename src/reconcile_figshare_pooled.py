#!/usr/bin/env python3
"""Reconcile the explicitly superseded mixed-density drafts to pooled-only drafts."""
import csv
import json
import sqlite3
from pathlib import Path

from upload_figshare_plots import Client


def main():
    root = Path(__file__).resolve().parents[1]
    old = root / "int/figshare_upload"
    new = root / "int/figshare_pooled_upload"
    old_links = list(csv.DictReader((old / "draft_links.csv").open()))
    records = list(csv.DictReader((new / "records.csv").open()))
    old_records = {r["record"]: r for r in csv.DictReader((old / "records.csv").open())}
    client = Client((root / "int/FS_TOKEN").read_text().strip())
    plan = []
    keep = set()
    for record in records:
        candidates = [r for r in old_links if
                      (r["record"] == record["record"] or r["record"].startswith(record["record"] + "__"))
                      and "Mean_cell_type_density" not in r["record"]]
        candidates.sort(key=lambda r: ("Pooled_synapse_density" not in r["record"], r["record"]))
        if not candidates:
            raise RuntimeError("No reusable draft for " + record["record"])
        chosen = candidates[0]
        if chosen["article_id"] in keep:
            raise RuntimeError("Duplicate reuse mapping")
        keep.add(chosen["article_id"])
        plan.append((record, chosen))
    # Preflight all 42 before mutation; only the known unpublished drafts qualify.
    details = {}
    for link in old_links:
        info = client.request("GET", "account/articles/" + link["article_id"])
        expected_title = "Synapse scatter plots — " + old_records[link["record"]]["title"]
        if info["title"] != expected_title or info.get("published_date") or info.get("is_public"):
            raise RuntimeError("Draft identity/publication check failed: " + link["article_id"])
        if link["article_id"] in keep and info["files"]:
            raise RuntimeError("Reusable draft unexpectedly contains files; inspect before regrouping")
        details[link["article_id"]] = dict(id=info["id"], title=info["title"], files=len(info["files"]))
    (new / "reconciliation_audit.json").write_text(json.dumps({
        "keep": sorted(keep), "original_drafts": details,
        "delete": [r["article_id"] for r in old_links if r["article_id"] not in keep],
    }, indent=2) + "\n")
    db = sqlite3.connect(new / "upload.sqlite")
    db.execute("CREATE TABLE IF NOT EXISTS records (name TEXT PRIMARY KEY, id INTEGER NOT NULL)")
    for record, chosen in plan:
        article_id = int(chosen["article_id"])
        client.request("PUT", f"account/articles/{article_id}", json={
            "title": "Synapse scatter plots — " + record["title"],
            "description": "<p>Individual synapse scatter and pooled-synapse depth-density plots. "
                "Only the original asis density mode is included. Filenames identify gene, neuropil, "
                "left/right hemisphere, presynaptic/postsynaptic compartment, and Notch condition. "
                "Legends are included. PDFs retain their original contents and are not merged.</p>",
        })
        db.execute("INSERT OR REPLACE INTO records VALUES (?, ?)", (record["record"], article_id))
        db.commit()
        print(f"Retained {article_id}: {record['title']}", flush=True)
    (old / "SUPERSEDED").write_text("Use ../figshare_pooled_upload. Mixed-density upload is retired.\n")
    for link in old_links:
        if link["article_id"] not in keep:
            client.request("DELETE", "account/articles/" + link["article_id"])
            with (new / "deleted_drafts.log").open("a") as handle:
                handle.write(link["article_id"] + "\n")
            print("Deleted unused draft " + link["article_id"], flush=True)
    print(f"Reconciliation complete: {len(keep)} retained, {len(old_links)-len(keep)} deleted", flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(str(exc) if isinstance(exc, RuntimeError) else type(exc).__name__, flush=True)
        raise SystemExit(1)
