#!/usr/bin/env python3
"""
Build a full seqkit-compatible rename map for an NCBI-flagged RagTag FASTA.

- The 32 chromosome/X/MT-assigned sequences (from known_map.tsv) get their
  fixed generic names (scaffold_01 ... scaffold_29, scaffold_X1, scaffold_X2, scaffold_MT).
- Every other header matching ^(NC|NW)_[0-9]+\.[0-9]+_RagTag$ (i.e. carries an
  embedded NCBI accession via RagTag's renaming) gets a generic sequential
  name: unplaced_scaffold_0001, unplaced_scaffold_0002, ...
- Everything else (ptg###### contigs -- hifiasm's own IDs, no NCBI accession
  embedded) is left unchanged and NOT written to the map (seqkit leaves
  unmapped headers alone).
- Anything that matches neither pattern is flagged to stderr for manual review
  rather than silently passed through.

Usage:
    python3 build_header_map.py <headers.txt> <known_map.tsv> <out_full_map.tsv>

headers.txt   : one header ID per line (first whitespace-delimited token, no '>')
known_map.tsv : the 32-row fixed map (old_id<TAB>new_id)
out_full_map.tsv : output, ready for `seqkit replace -k`
"""
import sys
import re

def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)

    headers_path, known_map_path, out_path = sys.argv[1:4]

    known_map = {}
    with open(known_map_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            old, new = line.split("\t")
            known_map[old] = new

    ragtag_accession_re = re.compile(r"^(NC|NW)_\d+\.\d+_RagTag$")

    headers = []
    with open(headers_path) as f:
        for line in f:
            h = line.strip()
            if h:
                headers.append(h)

    unplaced = []
    unchanged_ptg = 0
    flagged = []

    for h in headers:
        if h in known_map:
            continue  # handled directly from known_map
        elif ragtag_accession_re.match(h):
            unplaced.append(h)
        elif h.startswith("ptg"):
            unchanged_ptg += 1
        else:
            flagged.append(h)

    unplaced.sort()  # deterministic, reproducible numbering
    width = max(4, len(str(len(unplaced))))

    full_map = dict(known_map)  # start with the 32 fixed entries
    for i, h in enumerate(unplaced, start=1):
        full_map[h] = f"unplaced_scaffold_{i:0{width}d}"

    # sanity: every known_map old_id should actually be present in this file
    missing_known = [old for old in known_map if old not in headers]

    with open(out_path, "w") as out:
        for old, new in full_map.items():
            out.write(f"{old}\t{new}\n")

    # collision check
    new_names = list(full_map.values())
    dup_names = {n for n in new_names if new_names.count(n) > 1}

    print(f"Total headers in file:            {len(headers)}", file=sys.stderr)
    print(f"Chromosome/X/MT (known_map):       {len(known_map) - len(missing_known)} matched, {len(missing_known)} missing from file", file=sys.stderr)
    print(f"Unplaced RagTag-accession headers: {len(unplaced)} -> renamed", file=sys.stderr)
    print(f"ptg contigs left unchanged:        {unchanged_ptg}", file=sys.stderr)
    print(f"Unrecognized headers (NEEDS REVIEW): {len(flagged)}", file=sys.stderr)
    if flagged:
        print("  First 20 flagged headers:", file=sys.stderr)
        for h in flagged[:20]:
            print(f"    {h}", file=sys.stderr)
    if missing_known:
        print("  Known-map IDs not found in this file:", file=sys.stderr)
        for h in missing_known:
            print(f"    {h}", file=sys.stderr)
    if dup_names:
        print(f"  WARNING: {len(dup_names)} duplicate output names generated!", file=sys.stderr)

    print(f"\nWrote {len(full_map)} rename rules to {out_path}", file=sys.stderr)

if __name__ == "__main__":
    main()
