#author: haim krupkin+claude
#date: 11.13.2026
#description: step 10b. Pick example variants per cell type (3 each) for figures, by rules, from the main (20-PC) results.
# Then plot them: Rscript 10_plot_variant.R --file output/figures_epc20/variants/example_variants.csv
# Outputs: output/figures_epc20/variants/example_variants.csv, example_candidates_all_hits.csv (every hit with the criteria)
# rules: not a lineage-conflict gene; replicates pre/post (same sign as main); no other cell type with |z| > 2;
# the gene's expression rises with the cell type's proportion (Spearman rho >= 0.3), so the cell type really makes it
suppressMessages({library(dplyr); library(data.table)})
setwd("/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/rerun_LM22_aaron_1103/final_pipeline/output")
CT <- c("Neutrophils","Monocytes","Bcells","CD4","CD8","NK")
hits <- read.csv("figures_epc20/09_3_all_hits.csv") %>% filter(!lineage_conflict)
W_df <- read.csv("05_proportions_6classes.csv")
tpm <- readRDS("01_tpm_noHb.rds")
rd <- function(tag) rbindlist(lapply(unique(hits$chr), function(ch) fread(sprintf("06_scan_%s_top_variants%s.csv", ch, tag))[set == "genotype"]))
top <- rd("_epc20"); pre <- rd("_pre_epc20"); post <- rd("_post_epc20")
all_cand <- bind_rows(lapply(CT, function(ct) {
  h <- hits %>% filter(celltype == ct)
  t <- top[test == ct & gene_name %in% h$gene_name]
  z <- as.matrix(t[, paste0("b_", CT), with = FALSE]) / as.matrix(t[, paste0("se_", CT), with = FALSE])
  colnames(z) <- CT
  t$z_ct <- z[, ct]
  t$max_other_z <- apply(abs(z[, setdiff(CT, ct), drop = FALSE]), 1, function(r) if (all(is.na(r))) 0 else max(r, na.rm = TRUE))
  t$n_other_tested <- rowSums(!is.na(z[, setdiff(CT, ct), drop = FALSE]))
  b <- paste0("b_", ct)
  rep <- merge(unique(pre[, c("gene_name","variant_id", b), with = FALSE]), unique(post[, c("gene_name","variant_id", b), with = FALSE]),
               by = c("gene_name","variant_id"), suffixes = c("_pre","_post"))
  t <- merge(t, rep, by = c("gene_name","variant_id"))
  t$replicates <- sign(t[[paste0(b, "_pre")]]) == sign(t[[paste0(b, "_post")]]) & sign(t[[paste0(b, "_pre")]]) == sign(t[[b]])
  t$rho_expr_vs_prop <- sapply(t$gene_name, function(g) cor(log2(tpm[g, W_df$SAMPID] + 1), W_df[[ct]], method = "spearman"))
  t %>% as.data.frame() %>% transmute(gene_name, celltype = ct, variant_id, z = round(z_ct, 1), max_other_z = round(max_other_z, 1),
    n_other_tested, replicates, rho_expr_vs_prop = round(rho_expr_vs_prop, 2))
}))
# how many hits pass the expression check, per cell type (a quality signal in itself)
print(all_cand %>% group_by(celltype) %>% summarise(hits = n(), expr_rises_with_celltype = mean(rho_expr_vs_prop >= 0.3), .groups = "drop"))
pick <- all_cand %>% filter(replicates, max_other_z < 2, rho_expr_vs_prop >= 0.3) %>%
  group_by(celltype) %>% arrange(desc(abs(z))) %>% slice_head(n = 3) %>% ungroup() %>%
  mutate(celltype = factor(celltype, levels = CT)) %>% arrange(celltype) %>% as.data.frame()
print(pick)
write.csv(pick, "figures_epc20/variants/example_variants.csv", row.names = FALSE)
write.csv(all_cand, "figures_epc20/variants/example_candidates_all_hits.csv", row.names = FALSE)
