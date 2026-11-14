# cell_type_eQTLanalysis

Cell-type eQTLs (cis and trans) from bulk whole-blood RNA-seq (GTEx v8, 670 donors), without per-sample cell-type
profiles. Each genotype is tested inside one joint model with the cell-type proportions (TCA, Rahmani et al. 2019): all six
cell-type effects (neutrophils, monocytes, B, CD4, CD8, NK) are estimated together. The key elements:

- genotype effects are weighted by each cell type's share of the gene's RNA;
- the way-of-death batch effect is modelled;
- standard errors are robust (HC3);
- calls are calibrated against shuffled-genotype scans.

- `00_MASTER.sh`: runs the steps (`STEPS="..."`, `CHRS="..."`); `submit_master.sbatch` runs it as a SLURM job.
- `METHODS.md`: methods in plain language. `FILES.md`: every script and output. `TODO.md`: open items.
- Steps 1-5c: inputs, batch removal, CIBERSORTx proportions, RNA shares, hidden-noise PCs.
- Step 6: the cis scan (engine: `06_cis_scan.R` + compiled kernel `06_fit_variants.cpp`).
- Steps 7-11: calibrated calls and validation (pre/post replication, OneK1K, Kim-Hellmuth 2020, spike-ins, ASE, figures).
- Step 12 (`12a`, `12`, `12b`): cell-type trans-eQTLs. Candidate variants are tested against far genes with the same model,
  calibrated by shuffles, and checked for proportion-QTL artefacts.

This repository holds code and documentation only. Data (expression, genotypes, proportions, results) stays on the cluster
and is excluded by `.gitignore`. Paths in the scripts point to the Stanford SCG cluster.

Requirements: R 4.3 with dplyr, tidyr, data.table, snpStats, TCA, Rcpp, RcppArmadillo; plink2; the CIBERSORTx container
(step 4; credentials are passed as a file and never printed).
