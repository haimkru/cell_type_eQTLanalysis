# final_pipeline: to do (started 11.12.2026)

## Running / to report
- [x] Hidden-noise PCs (step 5c `05c_expression_pcs.R`, step 6 `N_EPC=10/20`), chr21-22 (11.12): hits about double
      (20 PCs: Neut 52->104, Mono 4->10, B 12->24, CD4 0->10, CD8 8->15, NK 0->10); 0 false calls per shuffle;
      pre/post same sign 100%; OneK1K enrichment holds (CD8 47% vs 15%, NK 60% vs 9%); main hits kept (effects rho
      0.95-1.00, same sign 100%). -> proposed for the main run (genome-wide: needs the user's go)
- [x] Genome-wide with N_EPC=20 (user's go 11.12; done 11.13, outputs tagged _epc20, figures in output/figures_epc20/):
      eGenes Neut 2345, Mono 146, B 99, CD4 86, CD8 317, NK 183 (was 1191/65/54/25/195/72); 0 false per shuffle;
      pre/post same sign 98-100%; raw-props rho 0.99; OneK1K enrichment p 4e-5 to 5e-55; 07b same direction 234/240;
      spike-in null 0/20; Kim-Hellmuth overlap OR 10.2 (p 1e-249), 56% of theirs in ours; ASE 8-26% vs 0-3% chance
- [x] _epc20 is the main result (user 11.13): N_EPC=20 default in step 6, steps 7-11 read _epc20 by default, METHODS.md / FILES.md / master updated

## Validation against Kim-Hellmuth et al. 2020 (Science; GTEx cell-type interaction eQTLs; PDFs in for_evin_cibersorting_expression/)
- [x] Overlap (`07c_kimhellmuth_check.R`, 11.12; their files in external/kim_hellmuth_ieqtl/): 31% of our neutrophil
      eGenes are their ieGenes vs 5% of other genes (OR 8.4, p 1.6e-143); 32% of theirs are ours; same lead signal
      (r2 >= 0.8) for 202/342 shared genes; direction agrees 80% (56/70) where our neutrophil-vs-lymphoid contrast is
      p < 0.05 (51% where it is not). Our neutrophil proportion vs their xCell score: rho 0.64.
- [x] ASE validation (`11_ase_validation.R`, 11.12; GTEx v8 phASER): eGenes whose heterozygotes' allelic ratio follows
      the cell type's proportion in the direction of our effect (p < 0.05), with allelic imbalance, vs shuffled:
      Neut 29% vs 2%, CD8 12% vs 1%, NK 13% vs 2%, CD4 21% vs 5%, Mono 8% vs 5%, B 9% vs 4%.
      Note: ASE can only see an effect that CHANGES with the proportion (other cell types also express the gene);
      for a marker gene made by one cell type only, the allelic ratio is the same in everyone -> lower ceiling.
- [x] 07c, 10, 11 in 00_MASTER.sh and FILES.md; as figures 11 and 12 in 09 (11.13)
- [ ] (user 11.13: not important for now) MAF filter within the low and high half of each cell type's share (theirs: MAF >= 5% in both halves), test on chr21-22
- [ ] (low priority) eigenMT gene-level correction (our M_eff test on chr21-22 gave no gain)

- [x] Leakage check (`07d_leakage_check.R`, 11.13, user's request): between-cell-type effect rho -0.18 to 0.58 (leaky profiles 0.96-0.99); structure: T/NK pairs share most, neutrophils most with monocytes, myeloid-lymphoid least; fine 10-pair order vs OneK1K not reproduced (rank 0.05; OneK1K summary matrix compressed 0.85-0.98); 89-92% of hits are bulk eGenes (48% of tested genes), 269 genes new vs bulk; hit sets mostly distinct

## Figures
- [x] Per-variant plot like their Fig. 1C: `10_plot_variant.R` (expression vs cell-type proportion, one line per
      genotype; any gene x cell type x variant; `--top N` or `--file list.csv` for many). Outputs in
      output/figures/variants/
- [x] Example variants (`10b_pick_example_variants.R`, 11.13): 3 per cell type by rules (no lineage conflict; replicates pre/post; no other cell type |z| > 2; expression rises with the cell type, rho >= 0.3) -> output/figures_epc20/variants/variants_example_variants.pdf. Found on the way: the monocyte proportion is the weakest estimate (classic markers CD14 rho 0.23, LYZ 0.25, VCAN 0.16 with it) -> monocyte labels least certain
- [x] Step 10 in 00_MASTER.sh (header) and FILES.md; plots now in output/figures_epc20/variants/ (main result)

## Trans-eQTLs (user 11.14: "cell type specific transqtls")
- [x] Steps 12a / 12 / 12b written; engine identical to step 6 (1e-15); chr22 test: 6 calls, 3 clean (~noise level), top
      call SIRPAP1 = pseudogene artefact -> built-in artefact flags added
- [ ] Genome-wide (STEPS="12 12b" CHRS=all 22; ~3-4 h on 16 CPUs) -- needs the user's go
- [ ] Optional: download the Saha & Battle 2018 cross-mappability list (public) and pass it as CROSSMAP_FILE
- [ ] Optional: GWAS blood-trait candidate variants (CAND=file) instead of / in addition to cis leads

## Waiting for the user's go
- [ ] edQTLs with the same design (user 11.13: later)
- [ ] 20 shuffles genome-wide (lower floor on the empirical p; ~30-40 h)
- [ ] Lab batch (SMGEBTCH) genome-wide (chr21-22: no change -> probably not needed)
