#!/bin/bash
#SBATCH --job-name=final_pipeline
#SBATCH --partition=interactive
#SBATCH --account=default
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/rerun_LM22_aaron_1103/final_pipeline/logs/slurm_%j.log
#author: haim krupkin+claude
#date: 11.07.2026
#description: runs 00_MASTER.sh as a SLURM job (keeps running if you disconnect).
#   all steps:   sbatch submit_master.sbatch
#   some steps:  sbatch --export=ALL,STEPS="7 8 9" submit_master.sbatch
echo "[$(date)] job ${SLURM_JOB_ID} on $(hostname); steps: ${STEPS:-all}"
bash /oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/rerun_LM22_aaron_1103/final_pipeline/00_MASTER.sh
