#!/bin/bash
#author: haim krupkin+claude
#date: 11.07.2026 (11.08: also the raw mixture, for the robustness check)
#description: step 4. CIBERSORTx Fractions (cell-type proportions), two runs:
#   main : mixture output/03_tpm_noHb_batch_removed.tsv (linear TPM, Hb-free, way-of-death batch removed)
#          -> output/04_cibersortx/fractions.tsv            (used by the eQTL scan)
#   raw  : mixture output/01_tpm_noHb_all_genes.tsv      (linear TPM, Hb-free, batch NOT removed)
#          -> output/04_cibersortx_raw/fractions.tsv        (robustness check only: do the eGenes hold?)
#   signature : LM22 (22 immune cell types)
#   settings  : B-mode on (platform correction), quantile normalisation off, 100 permutations, all 670 samples
#   runner    : run_cibersortx.py (+ cibersortx_native.py), the pinned CIBERSORTx container; credentials are read
#               from private/cibersortx_auth.json and never printed
# A run whose output folder already exists is skipped (delete the folder to rerun it). ~16 minutes per run.
set -euo pipefail

PIPE_DIR=/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/rerun_LM22_aaron_1103/final_pipeline
LM22=/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/deconvolution_LM22_20261006/reference/LM22.tsv
AUTH=/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/rerun_LM22_aaron_1103/private/cibersortx_auth.json
PYTHON=/home/hkrupkin/envs/catca-edit/bin/python
THREADS=${SLURM_CPUS_PER_TASK:-8}

## ================= one CIBERSORTx run =================
run_fractions () {
  local mixture=$1
  local out=$2
  if [ -d "${out}" ]; then
    echo "${out} exists; skipped (delete it to rerun)"
    return 0
  fi
  echo "[$(date)] CIBERSORTx Fractions start: ${mixture}"
  "${PYTHON}" "${PIPE_DIR}/run_cibersortx.py" \
    --stage fractions \
    --mixture "${mixture}" \
    --signature "${LM22}" \
    --outdir "${out}" \
    --auth "${AUTH}" \
    --threads "${THREADS}"
  echo "[$(date)] CIBERSORTx Fractions done: ${out}/fractions.tsv"
}

## ================= main: batch-removed mixture =================
run_fractions "${PIPE_DIR}/output/03_tpm_noHb_batch_removed.tsv" "${PIPE_DIR}/output/04_cibersortx"

## ================= robustness check: raw mixture =================
run_fractions "${PIPE_DIR}/output/01_tpm_noHb_all_genes.tsv" "${PIPE_DIR}/output/04_cibersortx_raw"
