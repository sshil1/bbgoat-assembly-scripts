import re

scaf_lens = {}
with open("/data/gmbluser/ncbi_annotation_prep/D863F/gene_qc/scaffold_lengths.tsv") as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 2:
            scaf_lens[parts[0]] = int(parts[1])

problem_genes = {}
with open("/data/gmbluser/ncbi_annotation_prep/D863F/D863F_dryrun_v2.val") as f:
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

gff_file = "/data/gmbluser/ncbi_annotation_prep/D863F/D863F_annotation_NCBI_v2.gff3"

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
internal_list = []

for gid, (scaf, start, end) in results.items():
    slen = scaf_lens.get(scaf, 0)
    if slen == 0:
        print(f"WARNING: scaffold {scaf} (gene {gid}) not in scaffold_lengths.tsv!")
        continue
    near_edge = (start <= 1000) or (slen - end <= 1000)
    span = end - start
    frac = span / slen if slen > 0 else 0
    is_main_chrom = scaf.startswith("scaffold_")
    large_fraction = frac > 0.3

    if near_edge or large_fraction:
        edge_count += 1
    else:
        internal_list.append((scaf, gid, start, end, slen, frac, sorted(problem_genes[gid]), is_main_chrom))

print(f"Edge/gap or large-scaffold-fraction (mark partial-eligible): {edge_count}")
print(f"Remaining candidates (main chrom + tiny fraction = TRUE internal errors): {len(internal_list)}")

with open("/data/gmbluser/ncbi_annotation_prep/D863F/gene_qc/problem_genes_final.tsv", "w") as out:
    out.write("scaffold\tgene_id\tstart\tend\tscaffold_len\tfraction\tproblem_types\tis_main_chrom\n")
    for scaf, gid, start, end, slen, frac, ptypes, is_main in sorted(internal_list, key=lambda x: not x[7]):
        out.write(f"{scaf}\t{gid}\t{start}\t{end}\t{slen}\t{frac:.4f}\t{','.join(ptypes)}\t{is_main}\n")

print("\n=== All internal candidates (main chromosome first) ===")
for row in sorted(internal_list, key=lambda x: not x[7]):
    print(row)
