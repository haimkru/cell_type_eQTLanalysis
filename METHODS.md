# Cell-type eQTLs from GTEx whole blood: methods (plain language)

Written 11.12.2026 from the scripts in this folder (00_MASTER.sh runs them). Updated 11.13: the main result is the
model WITH 20 hidden-noise PCs (step 5c; outputs tagged `_epc20`, figures in `output/figures_epc20/`). Numbers are from
that genome-wide run (11.13); the run without them (11.11) is kept untagged for comparison.

**Data.** GTEx v8 whole-blood RNA-seq and genotypes, 670 donors.

**Step 1, expression (01_expression_noHb.R).** Gene-level TPM. The 10 haemoglobin genes (HBA1, HBA2, HBB, HBD, HBG1,
HBG2, HBM, HBQ1, HBZ, HBE1) were removed, because they are a large and variable share of blood RNA, and each sample was
rescaled to 1 million. Genes with the same symbol were summed. For the eQTL analysis: genes with TPM > 0.1 in >= 20% of
samples.

**Step 2, covariates (02_covariates.R).** Per sample: blood draw before vs after death (pre-mortem = negative ischemic
time), ventilator case (Hardy scale 0), ischemic time, RIN, collection centre, sex, age.

**Step 3, remove the way-of-death effect from expression (03_remove_batch.R).** How the donor died is the largest source
of variation in GTEx blood (blood draw explains 70% of expression PC1). Each gene was fitted on blood draw, ventilator,
ischemic time, RIN and centre (quasi-Poisson, log link), and its expression divided by the fitted fold change. After:
0.2% of PC1. Pre-mortem, ventilator and centre B1 are largely the same donors, so their effects are removed together. The
corrected data are therefore not for studying pre- vs post-mortem differences.

**Step 4, proportions (04_cibersortx_fractions.sh).** CIBERSORTx Fractions, LM22 signature, on the corrected expression:
no quantile normalisation, B-mode batch correction, 100 permutations. A second run on the uncorrected expression is a
robustness check only.

**Step 5, 6 classes (05_proportions_6classes.R).** Neutrophils; monocytes (+ macrophages); B cells (naive, memory,
plasma); CD4 (all CD4 subsets + Tregs); CD8; NK (resting + activated). Small types left out (gamma-delta T, dendritic,
mast, eosinophils); rescaled to sum to 1, floor 1e-4. Means: neutrophils 54%, monocytes 16%, NK 10%, CD8 8%, CD4 7%,
B 5%.

**Step 5b, each cell type's share of each gene's RNA (05b_theta.R).** A cell type that is 10% of the cells can make 80%
of a gene's RNA, or none. Each gene's expression level per cell type was estimated by regressing its TPM on the 6
proportions across donors (non-negative least squares, floor 1e-3 = not expressed). From that, each cell type's share of
the gene's RNA per sample (phi). A cell type is tested in a gene only if its average share is >= 10%; otherwise NA.

**Step 5c, hidden-noise factors (05c_expression_pcs.R).** Bulk expression carries hidden technical noise (library
preparation, RNA quality beyond RIN, unknown batches); bulk eQTL studies remove it with expression PCs or PEER factors.
Plain expression PCs in blood mostly capture cell composition, which is the variation the model learns cell-type effects
from, so composition is removed first: each gene's log2(TPM + 1) is regressed on the 6 proportions, proportions x way of
death, and the known covariates (step 6 design without genotype), and PCA is run on what is left (genes scaled). The top
20 PCs are used (PC1 explains 16% of that remainder; the PCs are uncorrelated with the proportions). Tested first on
chr21-22 (10 and 20 PCs), then genome-wide: about twice the eGenes, all checks below unchanged.

**Step 6, the eQTL model (06_cis_scan.R, kernel 06_fit_variants.cpp).**
- One joint model per gene and variant: a separate genotype effect for every tested cell type, estimated together (based
  on TCA, Rahmani et al. 2019). Outcome log2(TPM + 1). Each cell type's genotype effect is weighted by its RNA share
  (phi); cell-type effects are learned from how proportions vary between people.
- Variants within 100 kb of the transcription start site, minor allele frequency >= 5%.
- Covariates: way of death (blood draw, ventilator) with its own effect in each cell type; genotype x way of death (so a
  cell-type effect cannot come from a difference between death groups); ischemic time, RIN, centre, sex, age, 5 genotype
  principal components, 20 hidden-noise PCs (step 5c).
- Weighted least squares, each sample weighted by its expected variance from a TCA fit without genotype; robust (HC3)
  standard errors, because a few samples with extreme proportions otherwise have too much influence.
- Tests: one per cell type (any genotype effect?); heterogeneity (do the effects differ?); neutrophil vs mean lymphoid.
- Gene-level p: smallest p over the gene's variants x number of variants (Bonferroni).
- Null: the whole scan repeated 5 times with donors' genotypes shuffled (all variants of a donor move together, so LD is
  kept; the genotype-expression link is broken).

**Step 7, eGene calls and checks (07_results_summary.R).**
- Calibration: each gene's p compared with the pooled gene-level p of all genes in the 5 shuffled scans (empirical p =
  share of null values at least as small); then Benjamini-Hochberg across genes; eGene = FDR < 5%.
- Checks: (1) each shuffled scan scored as real data -> 0 eGenes in every cell type; (2) eGenes refitted at the same
  variant in pre- and post-mortem donors separately -> same direction 98-100%, Spearman rho 0.90-0.94; (3) refitted with
  proportions from the uncorrected run -> rho 0.99, same direction 100%; (4) OneK1K (single-cell eQTLs, 982 donors):
  our eGenes are more often OneK1K eGenes in the matching cell type than other tested genes, in every cell type (Fisher
  p 4e-5 to 5e-55; e.g. NK 49% vs 8%, CD8 39% vs 11%, CD4 59% vs 18%; neutrophils compared with classical monocytes,
  as OneK1K is PBMC only).

**Step 7b, same signal as OneK1K? (07b_onek1k_lead_check.R).** For eGenes also in OneK1K: LD (r2, GTEx genotypes)
between our lead variant and OneK1K's (hg19 -> hg38 liftOver), and direction of effect (alleles matched by allele
frequency; frequencies near 50% left out). 239/470 genes r2 >= 0.8; same direction in 234/240 shared signals.

**Step 7c, GTEx neutrophil interaction eQTLs (07c_kimhellmuth_check.R).** Kim-Hellmuth et al. 2020 (Science) mapped
neutrophil interaction eQTLs in the same GTEx blood samples with a different method (xCell neutrophil score x genotype,
one cell type at a time, PEER factors, eigenMT; public files in external/kim_hellmuth_ieqtl). 56% of their 1,073
neutrophil ieGenes are among our neutrophil eGenes (27% of ours are theirs, vs 4% of other genes; odds ratio 10.2,
p 1e-249); 347/603 genes called by both have the same lead signal (r2 >= 0.8); where our neutrophil-vs-lymphoid contrast
is p < 0.05, their interaction has the same direction in 71% (132/185). Our neutrophil proportion vs their xCell score:
Spearman rho 0.64.

**Step 7d, leakage (07d_leakage_check.R).** Per-sample cell-type profiles (CIBERSORTx HiRes, TCA) had failed here: their
effects correlated 0.96-0.99 between all cell types and with bulk. (A) For each pair of cell types, the effects of the
genes called in either (same variant) correlate -0.18 to 0.58, not ~1; neutrophils vs B, CD4 and NK are about 0. The 6
effects are estimated together from correlated proportions, so the estimation alone already correlates them (shuffled
genotypes: -0.2 to -0.8 for most pairs, but +0.40 for neutrophils-CD8, which explains most of its real 0.57); measured
against that baseline, the T/NK pairs share the most (CD8-NK, CD4-CD8, CD4-NK) and neutrophil-lymphoid pairs the least,
and neutrophils share most with monocytes. Not as expected: monocytes vs CD8 and NK share about as much as B vs T cells,
which fits the weak monocyte proportion (caveat 5); the subtraction of the shuffled baseline is only a rough guide (the
shuffled values are from different genes). The finer order of the 10 PBMC pairs was not reproduced against OneK1K (Spearman 0.05; from OneK1K's
lead-variant summary files its own pairs only span 0.85-0.98, so there is little to rank). (B) At the same lead signal as
GTEx bulk, each cell type's effect correlates 0.88-0.96 with the bulk slope; this is expected for real eQTLs (bulk is the
sum of the cells) and does not separate leakage from real effects. (C) 89-92% of each cell type's eGenes are GTEx bulk
eGenes (vs 48% of all tested genes; same samples), 269 genes are new vs bulk, and the cell types' eGene sets are mostly
distinct (e.g. 207 of 317 CD8 eGenes are called only in CD8).
(D, 11.14) r_b (Qi et al. 2018; genes called in cell type A, their lead variant, effects in A and B, corrected for the
estimation-error correlation re taken from the shuffled scans; 09_results_figures.R section 13 run of 11.13, outputs
figures_epc20/09_13*_rb_*.csv): T/NK pairs 0.59-0.92, B vs T/NK 0.69-0.94, neutrophils vs B/CD4/NK -0.05 to 0.23,
neutrophils vs monocytes 0.37 (monocytes vs neutrophils 0.75). Exception: neutrophils vs CD8 0.70 (CD8 vs neutrophils
0.68), the pair with by far the largest re (+0.59; all others -0.29 to +0.35; their proportions correlate -0.72). Check
(lead variants of eGenes at FDR < 0.05, genome-wide): the slope of the CD8 z on the neutrophil z is 0.20 (354 genes),
and of the neutrophil z on the CD8 z 0.48 (65 genes), both below re, i.e. no more than the shared estimation error alone
gives; 46 neutrophil eGenes have a same-sign CD8 |z| > 3 (0 opposite). So the neutrophil-CD8 r_b is inflated by an error
correlation the correction cannot fully remove, and an effect called in one of the two cannot be cleanly assigned
between them; the other pairs are unaffected.

**Step 8, power (08_spikein_plan.R, 08_spikein_eval.R).** On chr22, genes without signal got a known genotype effect in
one cell type (0.5, 1, 2 or 4 SD per allele; 10 genes per cell type and size; 20 genes with no effect). 0/20 no-effect
genes called; found effects always in the right cell type; 50-90% found at 2 SD, 90-100% at 4 SD.

**Step 9, figures (09_results_figures.R).** All of the above (output/figures_epc20/). For the model without hidden-noise
PCs (output/figures/) also a lab-batch check: adding the GTEx expression batch (SMGEBTCH, 19 groups) as a covariate on
chr21-22 left the effects unchanged (rho 0.94-1.00, all same direction).

**Step 10, example variants (10_plot_variant.R).** Any gene x cell type x variant, as Kim-Hellmuth Fig. 1C: expression
(covariates removed) against the cell type's proportion, one regression line per genotype. Runs in RStudio or Rscript.

**Step 11, allele-specific expression (11_ase_validation.R).** Independent of the deconvolution: within a person, the two
copies of a gene are in the same cells, so their ratio is free of composition and batch differences between people. In
heterozygotes at each eGene's lead variant (GTEx phASER haplotype counts, >= 8 reads, >= 10 heterozygotes, Hampel outliers
removed), the allelic fold change should follow the cell type's proportion, in the direction of our effect. Out of all
eGenes (those ASE cannot test count as not validated), validated (Pearson p < 0.05, same direction) vs chance (the same
eGenes, proportions shuffled): neutrophils 19% (449/2,345) vs 2%, NK 11% (20/183) vs 1%, CD8 9% (28/317) vs 2%,
monocytes 8% (11/146) vs 2%, CD4 7% (6/86) vs 0%, B 4% (4/99) vs 0%. This can only see effects that change
with the proportion; for a gene made by one cell type only, the allelic ratio is the same in everyone, so many real hits
cannot pass it.

**Step 12, cell-type trans-eQTLs (12a, 12, 12b; added 11.14).**

*Candidates.* Testing every variant against every gene is about 1,000 times the cis scan and has little power at
n = 670, so only candidate variants are tested. They are the lead variants of the calibrated cell-type cis-eGenes (3,032;
2,709 after LD pruning at r2 < 0.8 within 1 Mb). Any other list, such as GWAS blood-trait variants, can be given instead.

*Model.* Each candidate is tested against every gene on another chromosome or more than 5 Mb away, with exactly the
step-6 model: phi weights and phi gate, 20 hidden-noise PCs, HC3 errors, the same kernel, and step 6's per-gene noise fit.
At pairs both scans cover, the results match step 6 to about 1e-15.

*Calibration.* The null is 5 shuffled-genotype scans: whole rows are shuffled, which keeps the LD between candidates and
the correlation between genes. The empirical FDR at a threshold is the mean number of shuffled pairs below it divided by
the number of real pairs below it.

*Checks.* Each call is checked for the known trans artefacts:
- the variant is a proportion QTL: lm(proportion ~ genotype + covariates), Bonferroni over called variants;
- hub variants with 10 or more trans genes;
- the trans gene is a pseudogene;
- the trans gene's name contains a cis gene's name, or the reverse (paralogs, pseudogenes);
- immunoglobulin, TCR, KIR or HLA genes;
- the trans gene cross-maps with any gene within 1 Mb of the variant (Saha & Battle 2018; GRCh38 / GENCODE v26 list,
  figshare 7315136; their recommended filter).

*chr22 test.* 706 genes, 7 minutes on 16 CPUs, so genome-wide is about 3-4 hours.
- The real scan kept no more pairs at p < 1e-4 than the shuffles (409 vs 394-429).
- Calibration: 0.2-0.4 false calls per test per shuffle.
- 6 calls:
  - the strongest (SIRPAP1, p 5e-44) is a mapping artefact: a processed pseudogene of SIRPA, at SIRPA's own cis lead;
  - one is an immunoglobulin gene;
  - one is at a proportion-QTL variant;
  - XPNPEP3 (NK and het) cross-maps with a gene near its variant;
  - only 1 is clean (RP4-669P10.16 in B), which is at the noise level (8 tests × 0.2).
- So trans effects at this sample size are expected to be few, and the artefact filters are needed.

**Result.** Calibrated eGenes: neutrophils 2,345, monocytes 146, B 99, CD4 86, CD8 317, NK 183 (without hidden-noise PCs:
1,191, 65, 54, 25, 195, 72).

**Caveats.** (1) 9 of 3,176 hits are receptor genes called outside their lineage (6 antibody genes in CD8, 3 T-cell
receptor genes in B / NK), flagged in figures_epc20/09_3_all_hits.csv. (2) Five shuffles put a floor of ~1e-5 on the
empirical p: rank top hits by the model p. (3) The corrected proportions are for the eQTL model only. (4) A hidden-noise
PC could partly carry a genetic effect, and adjusting for it could create spurious associations that shuffles cannot
detect; the checks with real genotypes (replication, OneK1K, Kim-Hellmuth, ASE) did not show this. (5) The monocyte
proportion is the weakest of the six estimates: even classic monocyte markers follow it only loosely (Spearman rho with
the proportion: CD14 0.23, LYZ 0.25, VCAN 0.16, FCN1 0.53, CD68 0.42). Monocyte eGenes are real eQTLs (OneK1K monocytes
28% vs 4%), but their monocyte label is the least certain.
