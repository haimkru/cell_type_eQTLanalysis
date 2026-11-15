#!/bin/bash
#author: haim krupkin+claude
#date: 11.07.2026 (11.08: one engine; old steps 6-9 moved to superseded/; calibrated calls; robustness check; compiled kernel)
manual review started by hk on 11/14/2026 -  end of  reviewed: 11/14/2026 code is good :)
#we also check we arent publishing anything confidental.
#description: MASTER script. Cell-type eQTLs from GTEx whole-blood bulk RNA-seq, without leakage and without the
# way-of-death batch effect. Every script of the pipeline is in this folder; outputs go to output/, logs to logs/.
#
# WHY THIS DESIGN (short):
#   - Per-sample cell-type profiles (CIBERSORTx HiRes, TCA tensor) only copy the bulk: their eQTLs were the bulk eQTL
#     in every cell type (effect correlation 0.96-0.99). So no profiles here.
#   - Instead the genotype is tested inside ONE model with the proportions (TCA, step 6): each cell type gets its own
#     genotype effect, all 6 estimated together, learned from how proportions vary between people.
#   - The biggest batch effect is the way of death (blood drawn before vs after death, ventilator). It shifts the
#     proportions, so it is removed before CIBERSORTx (step 3) and modelled again in the eQTL test (step 6).
#   - Robust (HC3) standard errors: plain ones gave false hits for small cell types (checked with shuffled genotypes).
#   - eGene calls are calibrated against 5 shuffled-genotype scans (empirical null), so the FDR holds by construction.
#
# STEPS                                   SCALE OF THE DATA                 MAIN OUTPUT (in output/)
#   1  01_expression_noHb.R               linear TPM, Hb genes removed      01_tpm_noHb_all_genes.tsv
#   2  02_covariates.R                    -                                 02_covariates.csv
#   3  03_remove_batch.R                  linear TPM (log-link fit)         03_tpm_noHb_batch_removed.tsv
#   4  04_cibersortx_fractions.sh         linear TPM, LM22                  04_cibersortx/fractions.tsv (+ _raw: check)
#   5  05_proportions_6classes.R          proportions                       05_proportions_6classes.csv (+ _raw)
#   5b 05b_theta.R                        linear TPM + proportions          05b_theta.csv, 05b_mean_phi.csv
#        (each gene's expression per cell type -> phi = the cell type's share of the gene's RNA, used by step 6)
#   5c 05c_expression_pcs.R               log2(TPM + 1)                     05c_expression_pcs.csv
#        (hidden-noise PCs: PCs of the expression left after removing composition + known covariates)
#   6  06_cis_scan.R  (per chromosome)    log2(TPM + 1), Hb-free            06_scan_<chr>_{genes,top_variants}_epc20.csv
#        11.10: genotype effects weighted by RNA share (phi), a cell type tested only where mean phi >= 0.10;
#        cis window 100 kb. 11.13: MAIN MODEL with 20 hidden-noise PCs (N_EPC=20, outputs tagged _epc20; steps 7-11
#        read _epc20 by default). Without them: N_EPC=0 OUT_TAG="" SCAN_TAG="" (11.11 results, untagged)
#        genotypes cut per chromosome; main model on all samples + 5 shuffles; then at the main best variants:
#        pre- and post-mortem draws separately (replication) and raw-mixture proportions (robustness)
#   7  07_results_summary.R               -                                 07_egenes.csv, 07_calibration.csv,
#        calibrated eGenes; checks: calibration, replication, robustness,    07_replication.csv, 07_onek1k.csv,
#        OneK1K; cell-type-specific genes                                    07_robustness_raw_proportions.csv,
#                                                                            07_specific_genes.csv
#   7b 07b_onek1k_lead_check.R            -                                 07b_onek1k_lead_summary.csv
#        eGenes also in OneK1K: same lead signal (LD r2, GTEx genotypes) and same direction?
#   8  08_spikein_plan.R + 06 (MODEL=spike)  log2(TPM + 1), chr22         08_spikein_power.csv
#      + 08_spikein_eval.R                  power check: known cell-type effects added, found?  (not in default STEPS)
#   7c 07c_kimhellmuth_check.R            -                                 07c_kimhellmuth_summary_epc20.csv
#        neutrophil eGenes vs the GTEx neutrophil interaction eQTLs of Kim-Hellmuth et al. 2020
#   7d 07d_leakage_check.R                -                                 07d_leakage_*_epc20.csv, figures_epc20/07d_*
#        leakage: effect correlation between cell types (vs shuffled + OneK1K), effect vs bulk, UpSet with bulk
#   9  09_results_figures.R               -                                 figures_epc20/09_all_figures.pdf
#        all figures (needs steps 5-8 and 7b of all 22 chromosomes; not in default STEPS)
#   10 10_plot_variant.R (by hand, also in RStudio)                         figures_epc20/variants/
#        any gene x cell type x variant: expression vs cell proportion, one line per genotype
#   11 11_ase_validation.R                allele-specific expression        11_ase_validation_summary_epc20.csv
#        (GTEx phASER) does the allelic ratio of heterozygotes follow the cell type's proportion? (needs step 9)
#   12 12a_trans_candidates.R, 12_trans_scan.R (per chromosome), 12b_trans_summary.R          12b_trans_calls_epc20.csv,
#        (11.14) cell-type trans-eQTLs: cis-eGene lead variants x all genes > 5 Mb away / other   12b_trans_egenes_epc20.csv,
#        chromosomes, the step-6 model; calibrated by shuffles; proportion-QTL / hub checks       12b_trans_calibration_epc20.csv
#        (not in default STEPS: STEPS="12a 12 12b")
# Helpers: run_cibersortx.py + cibersortx_native.py (step 4), 06_fit_variants.cpp (step 6). Files: FILES.md
#
# USAGE
#   run everything as a SLURM job:   sbatch submit_master.sbatch
#   run some steps:                  STEPS="6 7" bash 00_MASTER.sh
#   chromosomes for steps 6-7:       CHRS="chr21 chr22" (default); genome-wide: CHRS="chr1 chr2 ... chr22"
#   genome-wide (11.13, phi + 100 kb + 20 PCs): about 13 h for step 6 on 12 CPUs -> 4 jobs of STEPS=6 (chr1-5, 6-11,
#    12-16, 17-22), then STEPS="7 7b 7c 8 9 11" with all 22 in CHRS; chain them with --dependency=afterany
#   step 4 skips a CIBERSORTx run whose output folder exists: delete it to rerun.
set -euo pipefail

## ================= settings =================
#these are just general settings and arent that important
PIPE_DIR=/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/rerun_LM22_aaron_1103/final_pipeline
LOGS=${PIPE_DIR}/logs
GENO=/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/genotypes/GTEx_v8_WholeBlood_670donors_MAF01_plink1
PLINK2=/home/hkrupkin/envs/catca-edit/bin/plink2
#this enables skipping of stuff and such :)
STEPS=${STEPS:-"1 2 3 4 5 5b 5c 6 7"}
export CHRS=${CHRS:-"chr21 chr22"}
#we should awlays set it to not 0 as this increases accuracy.
# 11.13: main result = 20 hidden-noise PCs (tag _epc20; also names the logs). Set them to "" (and N_EPC=0) for the
# model without; ${X-default} keeps an explicitly empty value
export N_EPC=${N_EPC-20}
export OUT_TAG=${OUT_TAG-_epc20}
export SCAN_TAG=${SCAN_TAG-_epc20}
export NULL_TAG=${NULL_TAG-_epc20}
export SPIKE_TAG=${SPIKE_TAG-_epc20}
#this is to load R for the r steps
module load R/4.3.3 2>/dev/null || true
mkdir -p "${PIPE_DIR}/output"
mkdir -p "${LOGS}"

## ================= helper: run one step, write its log =================
#this helps us run steps while logging them, incase thigns fail :) 
run_step () {
  local name=$1
  shift
  echo "[$(date)] start ${name}"
  "$@" > "${LOGS}/${name}${OUT_TAG:-}.log" 2>&1
  echo "[$(date)] done  ${name}  (log: ${LOGS}/${name}${OUT_TAG:-}.log)"
}

## ================= steps =================
#this clarifies what each step does.
for s in ${STEPS}
do
  case ${s} in
    1) run_step 01_expression_noHb Rscript "${PIPE_DIR}/01_expression_noHb.R" ;;
    2) run_step 02_covariates Rscript "${PIPE_DIR}/02_covariates.R" ;;
    3) run_step 03_remove_batch Rscript "${PIPE_DIR}/03_remove_batch.R" ;;
    4) run_step 04_cibersortx_fractions bash "${PIPE_DIR}/04_cibersortx_fractions.sh" ;;
    5) run_step 05_proportions_6classes Rscript "${PIPE_DIR}/05_proportions_6classes.R" ;;
    # 5b (11.10): each gene's expression per cell type (theta) -> RNA-share weights (phi) used by step 6
    5b) run_step 05b_theta Rscript "${PIPE_DIR}/05b_theta.R" ;;
    # 5c (11.12): hidden-noise PCs (expression left after composition + covariates); used by step 6 when N_EPC > 0
    5c) run_step 05c_expression_pcs Rscript "${PIPE_DIR}/05c_expression_pcs.R" ;;
    # 7c / 11 (11.12): Kim-Hellmuth 2020 neutrophil comparison; allele-specific expression validation
    7c) run_step 07c_kimhellmuth_check Rscript "${PIPE_DIR}/07c_kimhellmuth_check.R" ;;
    7d) run_step 07d_leakage_check Rscript "${PIPE_DIR}/07d_leakage_check.R" ;;
    11) run_step 11_ase_validation Rscript "${PIPE_DIR}/11_ase_validation.R" ;;
    6)
      mkdir -p "${PIPE_DIR}/output/genotypes"
      for CHR in ${CHRS}
      do
        # genotypes of this chromosome (670 donors, MAF >= 1%), cut from the genome-wide GTEx plink file
        if [ ! -f "${PIPE_DIR}/output/genotypes/${CHR}.bed" ]; then
          "${PLINK2}" --bfile "${GENO}" --chr "${CHR}" --output-chr chrM --make-bed \
            --out "${PIPE_DIR}/output/genotypes/${CHR}" --threads 4 --memory 16000 > /dev/null
        fi
        # the main scan first (the checks below refit its best variants)
        CHR=${CHR} MODEL=deathadj run_step 06_cis_scan_${CHR} Rscript "${PIPE_DIR}/06_cis_scan.R"
        # replication: pre- and post-mortem draws separately
        CHR=${CHR} MODEL=pre run_step 06_cis_scan_${CHR}_pre Rscript "${PIPE_DIR}/06_cis_scan.R"
        CHR=${CHR} MODEL=post run_step 06_cis_scan_${CHR}_post Rscript "${PIPE_DIR}/06_cis_scan.R"
        # robustness: proportions from the raw (not batch-removed) mixture
        CHR=${CHR} MODEL=rawprops run_step 06_cis_scan_${CHR}_rawprops Rscript "${PIPE_DIR}/06_cis_scan.R"
      done
      ;;
    7) run_step 07_results_summary Rscript "${PIPE_DIR}/07_results_summary.R" ;;
    # 7b (11.11): stricter OneK1K check (lead-variant LD and direction), over CHRS
    7b) run_step 07b_onek1k_lead_check Rscript "${PIPE_DIR}/07b_onek1k_lead_check.R" ;;
    # 9 (11.11): all figures (always reads all 22 chromosomes)
    9) run_step 09_results_figures Rscript "${PIPE_DIR}/09_results_figures.R" ;;
    8)
      # power check (spike-in), on chr22 only: known cell-type effects added to real expression, then step 6 and scoring
      # (needs the main step-6 scans of all chromosomes: their shuffles are the null)
      CHR=chr22 run_step 08_spikein_plan Rscript "${PIPE_DIR}/08_spikein_plan.R"
      CHR=chr22 MODEL=spike run_step 06_cis_scan_chr22_spike Rscript "${PIPE_DIR}/06_cis_scan.R"
      CHR=chr22 run_step 08_spikein_eval Rscript "${PIPE_DIR}/08_spikein_eval.R"
      ;;
    # 12a / 12 / 12b (11.14): cell-type TRANS-eQTLs. 12a: candidate variants (lead variants of the calibrated cell-type
    # cis-eGenes, all 22 chromosomes, LD-pruned); 12: each candidate x every far gene of each chromosome in CHRS (same
    # model as step 6); 12b: calibrated calls + artefact checks over CHRS
    12a) CHRS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22" \
           run_step 12a_trans_candidates Rscript "${PIPE_DIR}/12a_trans_candidates.R" ;;
    12)
      for CHR in ${CHRS}
      do
        CHR=${CHR} run_step 12_trans_scan_${CHR} Rscript "${PIPE_DIR}/12_trans_scan.R"
      done
      ;;
    12b) run_step 12b_trans_summary Rscript "${PIPE_DIR}/12b_trans_summary.R" ;;
    *)
      echo "unknown step ${s}"
      exit 1
      ;;
  esac
done
echo "[$(date)] all requested steps done: ${STEPS}"
