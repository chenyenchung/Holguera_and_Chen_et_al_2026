#!/usr/bin/env python3
"""Stage unchanged CAM/TF PDFs under readable names for Figshare upload."""

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


NEUROPILS = {"ME": "Medulla", "LO": "Lobula", "LOP": "Lobula_plate"}
SIDES = {"L": "Left", "R": "Right"}
SYNAPSES = {"pre": "Presynaptic", "post": "Postsynaptic"}
DENSITIES = {"asis": "Pooled_synapse_density", "pertype": "Mean_cell_type_density"}
PATTERN = re.compile(
    r"^(ME|LO|LOP)_([LR])_(pre|post)_(asis|pertype)_(.+)\.pdf$"
)


def collect(root, density_mode="asis"):
    rows = []
    for family in ("CAM", "TF"):
        directory = root / "int" / f"P15_{family}"
        files = sorted(directory.rglob("*.pdf"))
        if not files:
            raise ValueError(f"No PDFs found in {directory}")
        for source in files:
            match = PATTERN.fullmatch(source.name)
            if not match:
                raise ValueError(f"Unrecognized filename: {source}")
            neuropil, side, synapse, density, tail = match.groups()
            if density != density_mode:
                continue
            if tail == "legend":
                gene, notch, kind = "", "", "legend"
            else:
                gene, notch = tail.rsplit("_", 1)
                if notch not in ("NotchOn", "NotchOff"):
                    raise ValueError(f"Unrecognized Notch condition: {source}")
                kind = "plot"
            prefix = [family, NEUROPILS[neuropil], SIDES[side], SYNAPSES[synapse]]
            suffix = ([gene, notch] if kind == "plot" else ["Legend"])
            filename = "__".join(prefix + suffix + [DENSITIES[density]]) + ".pdf"
            rows.append(dict(
                family=family, neuropil=NEUROPILS[neuropil],
                hemisphere=SIDES[side], synapse=SYNAPSES[synapse],
                density=DENSITIES[density], gene=gene, notch=notch,
                kind=kind, filename=filename,
                source=str(source.relative_to(root)), bytes=source.stat().st_size,
            ))
    return rows


def group_records(rows, limit):
    groups = defaultdict(list)
    for row in rows:
        groups[(row["family"], row["neuropil"], row["hemisphere"])].append(row)
    records = {}

    def split(key, entries, fields):
        if len(entries) <= limit:
            records["__".join(key)] = entries
        elif fields:
            subgroups = defaultdict(list)
            for row in entries:
                subgroups[row[fields[0]]].append(row)
            for value, subset in sorted(subgroups.items()):
                split((*key, value), subset, fields[1:])
        else:
            raise ValueError(f"Record {key} still exceeds {limit} files")

    for key, entries in sorted(groups.items()):
        split(key, entries, ["synapse", "density"])
    return records


def write_csv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path, default=Path("int/figshare_pooled_upload"))
    parser.add_argument("--max-files", type=int, default=500)
    args = parser.parse_args()
    root = args.root.resolve()
    output = args.output if args.output.is_absolute() else root / args.output
    if args.max_files < 1:
        parser.error("--max-files must be positive")
    if output.exists():
        parser.error(f"Output already exists; choose a new --output: {output}")
    rows = collect(root)
    records = group_records(rows, args.max_files)
    # Validate names and inputs before creating any staging files.
    names = [row["filename"] for row in rows]
    if len(names) != len(set(names)):
        raise ValueError("Renaming would produce duplicate filenames")
    for row in rows:
        with (root / row["source"]).open("rb") as handle:
            if handle.read(5) != b"%PDF-":
                raise ValueError(f"Not a PDF: {row['source']}")
    output.mkdir(parents=True)
    summaries = []
    manifest = []
    for record, entries in sorted(records.items()):
        folder = output / "records" / record
        folder.mkdir(parents=True)
        for row in sorted(entries, key=lambda r: r["filename"]):
            target = folder / row["filename"]
            target.symlink_to(root / row["source"])
            manifest.append(dict(record=record, staged_path=str(target.relative_to(output)), **row))
        summaries.append(dict(
            record=record, title=record.replace("__", " — ").replace("_", " "),
            files=len(entries), plots=sum(r["kind"] == "plot" for r in entries),
            legends=sum(r["kind"] == "legend" for r in entries),
            bytes=sum(r["bytes"] for r in entries),
        ))
    write_csv(output / "manifest.csv", manifest, list(manifest[0]))
    write_csv(output / "records.csv", summaries, list(summaries[0]))
    readme = f"""# Figshare upload staging

{len(rows):,} unchanged PDFs in {len(records)} proposed Figshare records.
Maximum PDFs per record: {max(r['files'] for r in summaries)}.

Nothing has been uploaded or published. Original filenames and PDF contents
are unchanged. Files under records/ are absolute symbolic links to the originals,
using readable upload names; they avoid duplicating the full dataset locally.
An uploader must open each link as a file and explicitly use the filename column
as the remote filename. To move staging to another machine, dereference links.

manifest.csv maps every source to its staged name, record, and plot conditions.
records.csv lists proposed record titles, file counts, and sizes. These two
indexes and this README are local planning files, not additional files in each
record. Records may already contain exactly 500 PDFs.

## Naming and organization

Filename fields are separated by double underscores:
family, neuropil, hemisphere, synapse compartment, gene, Notch condition,
density method. Legend files have Legend in place of gene and Notch condition.
Gene symbols are preserved verbatim. Spaces within descriptive fields use
underscores; Left and Right refer to anatomical hemisphere.

Only Pooled_synapse_density (original asis mode) is included: a density estimate
using pooled synapse depths. The pertype version is excluded.
No PDFs are merged. Legends are retained with their condition.

Records begin with family × neuropil × hemisphere. When over the file limit,
they split by presynaptic/postsynaptic, then by density method if still needed.
No alphabetical or gene-based record splits are needed with the current inputs.

## Upload prerequisites

Use the destination account's API token through a local secret/environment
variable, never in this manifest or source control. Account type and available
quota must be established before uploading. Draft metadata needs the intended
authors, description, categories, and license before publication. The script
only prepares local files and does not authenticate, upload, or publish.
"""
    (output / "README.md").write_text(readme)
    print(f"Prepared {len(rows):,} PDFs in {len(records)} records at {output}")
    print(f"Largest record: {max(r['files'] for r in summaries)} PDFs")


if __name__ == "__main__":
    main()
