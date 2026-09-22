#!/usr/bin/env bash
# Reconstructed from shell history (2026-08-31). Produces genes_to_exclude.txt
# from the categorize_problem_genes*.py output (problem_genes_v3.tsv / problem_genes_final.tsv).
set -euo pipefail

# --- 803M pipeline: adds a borderline-vs-true-internal distinction ---
awk -F'\t' 'NR==1{print; next} {dist_start=$3; dist_end=$5-$4; if(dist_start<=1000 || dist_end<=1000) print $0"\tBORDERLINE"; else print $0"\tTRUE_INTERNAL"}' \
    803M_problem_genes_v3.tsv > 803M_problem_genes_v3_widened.tsv
awk -F'\t' '$1 ~ /^scaffold_/ && $0 ~ /TRUE_INTERNAL/ {print $2}' \
    803M_problem_genes_v3_widened.tsv > 803M_genes_to_exclude.txt

# --- D863F pipeline: no borderline-widening step was applied ---
awk -F'\t' 'NR>1 {print $2}' \
    D863F_problem_genes_final.tsv > D863F_genes_to_exclude.txt
