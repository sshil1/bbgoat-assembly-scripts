# BBGOAT genome assembly — custom scripts

Scripts used in the Black Bengal Goat de novo genome assembly and annotation
(BioProject PRJNA1502507). See the associated Scientific Data manuscript for context.

- `build_header_map.py` — FASTA header renaming for NCBI submission
- `categorize_problem_genes_v3.py` — gene-model edge/length QC categorization for 803M
  (final version; supersedes earlier iterations not included here)
- `categorize_problem_genes_d863f.py` — gene-model edge/length QC categorization for D863F
- `gene_qc_exclude_list.sh` — post-processing step producing each individual's final
  `genes_to_exclude.txt`. Note: the 803M and D863F pipelines are not identical — 803M's
  list additionally excludes only genes classified TRUE_INTERNAL (>1000 bp from a
  scaffold edge) after a borderline-widening pass; D863F's exclude list was taken
  directly from the categorization script's output with no borderline-widening step.
