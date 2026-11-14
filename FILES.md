# final_pipeline: files

Folder: `/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/rerun_LM22_aaron_1103/final_pipeline/`

One pipeline, steps 1-11 (with 5b, 5c, 7b, 7c), run by `00_MASTER.sh`. Methods in plain language: `METHODS.md`; open items: `TODO.md`. Old versions of steps are in `superseded/` (not used).

**Main result (11.13): the model with 20 hidden-noise PCs; its outputs carry the tag `_epc20`** (e.g. `06_scan_chr1_genes_epc20.csv`, `07_egenes_epc20.csv`, `figures_epc20/`). Untagged files = the same model without those PCs (11.11), kept for comparison.

## Scripts

| File | What it does |
|---|---|
| `00_MASTER.sh` | runs the steps in order; `STEPS="..."` picks steps, `CHRS="..."` picks chromosomes |
| `submit_master.sbatch` | runs the master as a SLURM job (interactive partition, 16 CPUs) |
| `01_expression_noHb.R` | GTEx whole-blood TPM, haemoglobin genes removed, re-TPM, summed by symbol, expressed-gene subset |
| `02_covariates.R` | per-sample table: blood draw, ventilator, ischemic time, RIN, centre, sex, age |
| `03_remove_batch.R` | removes the way-of-death batch from the TPM (per-gene log-link fit, divide by the batch fold change) |
| `04_cibersortx_fractions.sh` | CIBERSORTx Fractions, LM22: main run on the batch-removed TPM; second run on the raw TPM (robustness check only) |
| `05c_expression_pcs.R` | hidden-noise PCs: PCs of the expression left after removing composition + known covariates (step 6 uses the top 20) |
| `05_proportions_6classes.R` | 22 LM22 types -> 6 classes (Neutrophils, Monocytes, B, CD4, CD8, NK), sum to 1; both runs |
| `05b_theta.R` | each gene's expression per cell type (theta, NNLS of linear TPM on the proportions) -> phi = the cell type's share of the gene's RNA (used by step 6) |
| `06_cis_scan.R` | THE eQTL engine: one joint TCA model per gene x variant, all 6 cell-type effects + contrasts, robust (HC3) errors, 5 shuffled-genotype null scans. Genotype effects weighted by RNA share (phi; a cell type is tested only where its mean phi >= 0.10, else NA); cis window 100 kb. MODEL = deathadj (main) / pre / post (replication) / rawprops (robustness) / spike (step 8) |
| `06b_meff.R` | (tested 11.10, not used) effective number of independent variants per gene (Li & Ji); gave no gain over Bonferroni |
| `06_fit_variants.cpp` | compiled kernel of step 6 (the per-variant regression; same math as the R version, ~4x faster) |
| `07_results_summary.R` | calibrated eGenes (empirical null from the shuffles); checks: calibration, pre/post replication, OneK1K, raw-proportions robustness; cell-type-specific genes |
| `07b_onek1k_lead_check.R` | stricter OneK1K check: for eGenes also OneK1K eGenes, LD (GTEx genotypes) between our lead and OneK1K's lead (hg19 -> hg38 liftOver) and the direction of effect |
| `07c_kimhellmuth_check.R` | our neutrophil eGenes vs the GTEx neutrophil interaction eQTLs of Kim-Hellmuth et al. 2020 (overlap, same lead signal, direction) |
| `07d_leakage_check.R` | leakage check: effect correlation between cell types (vs shuffled baseline and vs OneK1K), effect vs bulk, UpSet with GTEx bulk eGenes |
| `08_spikein_plan.R`, `08_spikein_eval.R` | step 8 (power check, chr22): known cell-type effects added to real expression (06 with MODEL=spike), then scored: detected? right cell type? leak? |
| `09_results_figures.R` | ONE script with all figures: proportions, batch effects, hits, real vs shuffled, replication, robustness, OneK1K overlap, spike-in power, OneK1K lead variants, lab batch (no-PC run only), Kim-Hellmuth (7c), ASE (11) |
| `10_plot_variant.R` | plot any gene x cell type x variant (expression vs cell proportion, one line per genotype); RStudio or Rscript (`--top N`, `--file list.csv`) |
| `10b_pick_example_variants.R` | picks 3 example variants per cell type by rules (no lineage conflict, replicates, specific, expression follows the cell type) for step 10 |
| `11_ase_validation.R` | allele-specific expression check (GTEx phASER): does the allelic ratio of heterozygotes follow the cell type's proportion? |
| `12a_trans_candidates.R` | (11.14) trans candidates: lead variants of the calibrated cell-type cis-eGenes (all chromosomes), LD-pruned (r2 < 0.8 within 1 Mb); or any variant list (`CAND=file`) |
| `12_trans_scan.R` | (11.14) cell-type trans-eQTL scan per chromosome of genes: every candidate x every gene on another chromosome or > 5 Mb away; the step-6 model and kernel, step-6 null fits reused, 5 shuffles; keeps pairs with p < 1e-4 + test counts |
| `12b_trans_summary.R` | (11.14) calibrated trans calls (empirical FDR vs shuffles), calibration, proportion-QTL and hub checks, cell-type specificity, cis link, optional cross-mappability filter |
| `run_cibersortx.py`, `cibersortx_native.py` | CIBERSORTx runner (pinned container; credentials never printed) |

## Outside inputs (read only)

| File | Used in |
|---|---|
| `expression/gtex_wholeblood_gene_tpm.csv.gz` | 1 |
| GTEx v8 `SampleAttributesDS`, `SubjectPhenotypesDS` (`catca-edit-hpc/edQTLs_GTEx/sample_annotations/`) | 2 |
| `deconvolution_LM22_20261006/reference/LM22.tsv` | 4 |
| `private/cibersortx_auth.json` (owner-only) | 4 |
| `genotypes/GTEx_v8_WholeBlood_670donors_MAF01_plink1.{bed,bim,fam}` (670 donors, MAF >= 1%, genome-wide) | 6 (cut per chromosome by the master) |
| GENCODE v26 GTF, GTEx v8 Whole_Blood eQTL covariates (genotype PCs 1-5) | 6 |
| `external/onek1k_top_eqtl/` (OneK1K top eQTL per gene, Zenodo 18910121) | 7 |
| `external/kim_hellmuth_ieqtl/` (public GTEx portal: `GTEx_Analysis_v8_ieQTL.tar`, xCell scores; downloaded 11.12) | 7c |
| GTEx v8 phASER (`/labs/smontgom/grps_smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/ase/GTEx_Analysis_v8_phASER`) | 11 |

## Outputs (`output/`)

| File | Step | Content |
|---|---|---|
| `01_tpm_noHb_all_genes.tsv`, `01_tpm_noHb.rds` | 1 | linear TPM, Hb-free, all genes |
| `01_gene_subset.txt` | 1 | expressed genes (TPM > 0.1 in >= 20% of samples) |
| `02_covariates.csv` | 2 | covariates per sample |
| `03_tpm_noHb_batch_removed.tsv` / `.rds` | 3 | batch-removed linear TPM (CIBERSORTx input) |
| `03_batch_check.csv` | 3 | PC1 ~ blood draw, before vs after |
| `04_cibersortx/fractions.tsv` | 4 | 22 LM22 proportions + fit per sample (main) |
| `04_cibersortx_raw/fractions.tsv` | 4 | the same on the raw TPM (robustness check) |
| `05_proportions_6classes.csv` | 5 | 6-class proportions (eQTL input). Way of death was removed from the expression before CIBERSORTx, so these are NOT for studying pre- vs post-mortem composition (use the raw run for that) |
| `05_proportions_6classes_raw.csv` | 5 | 6-class proportions from the raw run (robustness check) |
| `05_proportions_by_draw.csv` | 5 | mean per class, pre vs post-mortem, both runs |
| `genotypes/<chr>.{bed,bim,fam}` | 6 | genotypes of one chromosome |
| `06_scan_<chr>_genes{,_pre,_post,_rawprops}.csv` | 6 | per gene x test: best variant, number of variants, gene-level p, BH (real and shuffled) |
| `06_scan_<chr>_top_variants{,_pre,_post,_rawprops}.csv` | 6 | the same rows with all 6 effects, SEs and p-values |
| `06_scan_<chr>{,_pre,_post,_rawprops}_null_fit.rds` | 6 | TCA fit without genotype (sigma per gene and cell type, tau) |
| `07_egenes.csv` | 7 | calibrated eGenes per test (and the uncalibrated count for comparison) |
| `07_calibration.csv` | 7 | false eGenes when each shuffle is scored against the others |
| `07_replication.csv` | 7 | main eGenes in pre vs post-mortem draws: same sign, rank correlation |
| `07_robustness_raw_proportions.csv` | 7 | main eGenes refitted with the raw-run proportions |
| `07_onek1k.csv` | 7 | OneK1K overlap per cell type |
| `07_specific_genes.csv` | 7 | cell-type-specific eGenes (calibrated het test) and the cell type carrying the effect |
| `07_celltype_pattern.csv` | 7 | per cell type: are its eGenes' effects concentrated in that cell type? (nominal per-gene contrasts) |
| `08_spike_plan_chr22{,_with_variants}.csv` | 8 | spiked genes: cell type, effect size (SD per allele), variant |
| `06_scan_chr22_{genes,top_variants}_spike.csv` | 8 | step-6 scan of the spiked genes |
| `08_spikein_power.csv`, `08_spikein_genes.csv` | 8 | power per cell type x effect size; per-gene result |
| `05b_theta.csv`, `05b_mean_phi.csv` | 5b | expression per cell type per gene; mean RNA share per gene and cell type |
| `07b_onek1k_lead_pairs.csv`, `07b_onek1k_lead_summary.csv` | 7b | per eGene x OneK1K type: LD r2, direction; summary per cell type |
| `06_scan_chr2{1,2}_{genes,top_variants}_labbatch.csv`, `06_scan_chr2{1,2}_labbatch_null_fit.rds` | 6 (LABBATCH=1) | robustness test 11.11: + GTEx expression batch as covariate, chr21-22 only (results unchanged; figure 10) |
| `figures/09_<n>_<name>.{png,csv}`, `figures/09_all_figures.pdf`, `figures/09_3_all_hits.csv` | 9 | all figures, the data behind each, every calibrated hit |
| `*_epc20.*` (06, 07, 07b, 07c, 08, 11), `figures_epc20/` | 5c-11 | MAIN RESULT (11.13): the same outputs for the model with 20 hidden-noise PCs; `figures_epc20/variants/` = step-10 plots |
| `05c_expression_pcs.csv`, `05c_variance_explained.csv` | 5c | hidden-noise PCs per sample; variance each explains |
| `07c_kimhellmuth_{pairs,summary}*.csv` | 7c | per gene shared with Kim-Hellmuth: lead variants, LD, direction; summary |
| `11_ase_validation_{genes,summary}*.csv` | 11 | per eGene: heterozygotes, aFC vs proportion correlation; summary per cell type |
| `07d_leakage_{pairs,bulk,overlap,onek1k_structure}*.csv`, `figures_epc20/07d_*.png` | 7d | leakage check tables and figures |
| `superseded/window1Mb_1110/`, `superseded/phi_tests_chr21_22_1110/` | - | old genome-wide 1 Mb cell-share results; chr21-22 phi test summaries (kept, not used) |
