import re

# 1. Load scaffold lengths
scaf_lens = {}
with open("/data/gmbluser/ncbi_annotation_prep/803M/gene_qc/scaffold_lengths.tsv") as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 2:
            scaf_lens[parts[0]] = int(parts[1])

# 2. Extract broken gene IDs + problem types from .val file
problem_genes = {}
with open("/data/gmbluser/ncbi_annotation_prep/803M/803M_dryrun_v2.val") as f:
    for line in f:
        if "Error: valid [" not in line:
            continue
        ptype_m = re.search(r'\[(SEQ_[^\]]+)\]', line)
        gene_m = re.search(r'cds\.(g\d+)', line)
        if ptype_m and gene_m:
            gid = gene_m.group(1)
            pt = ptype_m.group(1)
            problem_genes.setdefault(gid, set()).add(pt)

print(f"Total unique problem genes found in .val file: {len(problem_genes)}")

# 3. Use the EXACT, correct GFF3 file (explicit, no fragile matching)
gff_file = "/data/gmbluser/ncbi_annotation_prep/803M/803M_annotation_NCBI_v2.gff3"

# 4. Parse GFF for exact coordinates (gene-level span, using 'gene' feature specifically)
results = {}
with open(gff_file) as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9 or parts[2] != "gene":
            continue
        scaf = parts[0]
        start = int(parts[3])
        end = int(parts[4])
        m_id = re.search(r'ID=(g\d+)', parts[8])
        if m_id:
            gid = m_id.group(1)
            if gid in problem_genes:
                results[gid] = (scaf, start, end)

print(f"Successfully mapped coordinates in GFF for: {len(results)}\n")

edge_count = 0
internal_count = 0
internal_list = []

for gid, (scaf, start, end) in results.items():
    slen = scaf_lens.get(scaf, 0)
    if slen == 0:
        print(f"WARNING: scaffold {scaf} (gene {gid}) not found in scaffold_lengths.tsv -- name mismatch!")
        continue
    near_edge = (start <= 100) or (slen - end <= 100)
    if near_edge:
        edge_count += 1
    else:
        internal_count += 1
        internal_list.append((scaf, gid, start, end, slen, sorted(problem_genes[gid])))

print(f"Near scaffold edge/gap (likely genuine partial): {edge_count}")
print(f"NOT near edge (likely real gene-calling error): {internal_count}")

with open("/data/gmbluser/ncbi_annotation_prep/803M/gene_qc/problem_genes_v3.tsv", "w") as out:
    out.write("scaffold\tgene_id\tstart\tend\tscaffold_len\tproblem_types\n")
    for scaf, gid, start, end, slen, ptypes in internal_list:
        out.write(f"{scaf}\t{gid}\t{start}\t{end}\t{slen}\t{','.join(ptypes)}\n")

print("\n=== Sample internal (real error) genes ===")
for row in internal_list[:15]:
    print(row)
