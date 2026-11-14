#author: haim krupkin+claude
#date: 11.13.2026
#description: step 7d. Leakage check: are the cell-type eQTLs distinct, or one bulk signal copied into every cell type?
# The CIBERSORTx HiRes / TCA profiles failed exactly here (effect correlation 0.96-0.99 between all cell types and with
# bulk; most hits in every cell type = bulk eGenes). Three checks, all on the main (20-PC) results:
#   A. do betas cross between cell types? per pair (A, B): genes called in A or in B (both cell types tested), at the
#      calling cell type's lead variant: Spearman rho of the A effect vs the B effect; and the share of A's eGenes that
#      also have a nominal effect in B (|z| > 2, same sign). Baseline: the same on the 5 shuffled scans (the 6 effects are
#      estimated together from correlated proportions, so the estimation alone gives some correlation)
#   B. effect vs bulk: each cell type's effect vs the GTEx whole-blood bulk eQTL slope, at the bulk lead variant when it
#      is the same signal as ours (LD r2 >= 0.8, GTEx genotypes; orientation exact, both ALT-coded)
#   C. discovery overlap: UpSet of the 6 cell types' eGenes + GTEx bulk whole-blood eGenes (qval <= 0.05)
#   D. OneK1K reference: the same between-cell-type effect correlation in OneK1K single cells; do the 10 pairs rank
#      the same (related cell types share more)? Leakage = all pairs alike and high
# Inputs : output/06_scan_<chr>_top_variants<SCAN_TAG>.csv, figures<SCAN_TAG>/09_3_all_hits.csv, output/genotypes/<chr>,
#          GTEx v8 Whole_Blood.v8.egenes.txt.gz
# Outputs: output/07d_leakage_pairs<SCAN_TAG>.csv (A), 07d_leakage_bulk<SCAN_TAG>.csv (B), 07d_leakage_overlap<SCAN_TAG>.csv
#          (C), figures<SCAN_TAG>/07d_leakage.png and 07d_upset_bulk.png
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(snpStats))

## ================= settings =================
pipe_dir <- "/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/rerun_LM22_aaron_1103/final_pipeline"
out_dir  <- file.path(pipe_dir, "output")
SCAN_TAG <- Sys.getenv("SCAN_TAG", "_epc20")    # main result = 20 hidden-noise PCs; "" = the model without
fig_dir  <- file.path(out_dir, paste0("figures", SCAN_TAG))
bulk_file <- "/labs/smontgom/grps_smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/eqtl/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.egenes.txt.gz"
CHRS <- paste0("chr", 1:22)
CELLTYPES <- c("Neutrophils", "Monocytes", "Bcells", "CD4", "CD8", "NK")
CT_COLS <- c(Neutrophils = "#1f77b4", Monocytes = "#ff7f0e", Bcells = "#2ca02c", CD4 = "#d62728", CD8 = "#9467bd", NK = "#8c564b")
theme_set(theme_bw(base_size = 12))

## ================= read =================
hits <- read.csv(file.path(fig_dir, "09_3_all_hits.csv"))
top <- bind_rows(lapply(CHRS, function(ch) fread(file.path(out_dir, paste0("06_scan_", ch, "_top_variants", SCAN_TAG, ".csv"))) %>%
  filter(test %in% CELLTYPES) %>%
  select(chr, gene_name, test, set, variant_id, p_gene, all_of(paste0("b_", CELLTYPES)), all_of(paste0("se_", CELLTYPES)))))

## ================= A. do betas cross between cell types? =================
# lead variant row of the calling cell type; real = the eGenes; null = each shuffle's 300 best genes of that cell type
lead_rows <- function(s, genes_by_ct) {
  bind_rows(lapply(CELLTYPES, function(ct) top %>% filter(set == s, test == ct, gene_name %in% genes_by_ct[[ct]])))
}
pair_stats <- function(rows, label) {
  bind_rows(lapply(CELLTYPES, function(a) bind_rows(lapply(setdiff(CELLTYPES, a), function(b) {
    d <- rows %>% filter(test %in% c(a, b)) %>% distinct(gene_name, .keep_all = TRUE)    # a gene called in both: once
    ba <- d[[paste0("b_", a)]]; bb <- d[[paste0("b_", b)]]
    za <- ba / d[[paste0("se_", a)]]; zb <- bb / d[[paste0("se_", b)]]
    ok <- !is.na(ba) & !is.na(bb)
    own <- rows %>% filter(test == a)                                                  # A's eGenes only
    zb_own <- own[[paste0("b_", b)]] / own[[paste0("se_", b)]]
    za_own <- own[[paste0("b_", a)]] / own[[paste0("se_", a)]]
    data.frame(data = label, celltype_a = a, celltype_b = b, genes_both_tested = sum(ok),
               rho = if (sum(ok) >= 10) cor(ba[ok], bb[ok], method = "spearman") else NA,
               share_a_egenes_nominal_in_b = mean(abs(zb_own) > 2 & sign(zb_own) == sign(za_own), na.rm = TRUE))
  }))))
}
genes_real <- lapply(setNames(CELLTYPES, CELLTYPES), function(ct) hits$gene_name[hits$celltype == ct])
pairs_real <- pair_stats(lead_rows("genotype", genes_real), "real genotypes")
null_sets <- sort(unique(top$set[top$set != "genotype"]))
pairs_null <- bind_rows(lapply(null_sets, function(s) {
  g <- lapply(setNames(CELLTYPES, CELLTYPES), function(ct)
    top %>% filter(set == s, test == ct) %>% arrange(p_gene) %>% slice_head(n = 300) %>% pull(gene_name))
  pair_stats(lead_rows(s, g), s)
})) %>%
  group_by(celltype_a, celltype_b) %>%
  summarise(data = "shuffled genotypes (mean of 5)", genes_both_tested = round(mean(genes_both_tested)),
            rho = mean(rho, na.rm = TRUE), share_a_egenes_nominal_in_b = mean(share_a_egenes_nominal_in_b, na.rm = TRUE),
            .groups = "drop")
pairs <- bind_rows(pairs_real, pairs_null) %>%
  mutate(celltype_a = factor(celltype_a, levels = CELLTYPES), celltype_b = factor(celltype_b, levels = CELLTYPES))
cat("A. between-cell-type effect correlation (rho, symmetric part), real vs shuffled:\n")
print(pairs %>% filter(as.integer(celltype_a) < as.integer(celltype_b)) %>%
        select(data, celltype_a, celltype_b, genes_both_tested, rho) %>%
        pivot_wider(names_from = data, values_from = c(genes_both_tested, rho)) %>% as.data.frame(), digits = 2)

pA <- ggplot(pairs %>% filter(as.integer(celltype_a) < as.integer(celltype_b)),
             aes(celltype_b, celltype_a, fill = rho)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = ifelse(is.na(rho), "", sprintf("%.2f", rho))), size = 3.5) +
  facet_wrap(~ data) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-1, 1), na.value = "grey90") +
  labs(x = NULL, y = NULL, fill = "rho",
       title = "A. Do effects cross between cell types? (genes called in either; same variant)",
       subtitle = "leaky CIBERSORTx/TCA profiles gave 0.96-0.99 for every pair; shuffled = what the joint estimation alone gives") +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1))
pA2 <- ggplot(pairs %>% filter(data == "real genotypes"), aes(celltype_b, celltype_a, fill = share_a_egenes_nominal_in_b)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = scales::percent(share_a_egenes_nominal_in_b, accuracy = 1)), size = 3.5) +
  scale_fill_gradient(low = "white", high = "#b2182b", limits = c(0, 1), labels = scales::percent) +
  labs(x = "also an effect in (|z| > 2, same sign)", y = "eGenes of", fill = NULL,
       title = "Share of a cell type's eGenes with an effect in another cell type") +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1))

## ================= B. effect vs bulk effect (same signal) =================
bulk <- fread(bulk_file, select = c("gene_name", "variant_id", "slope", "qval")) %>%
  group_by(gene_name) %>% filter(n() == 1) %>% ungroup() %>%
  rename(bulk_variant = variant_id, bulk_slope = slope, bulk_qval = qval)
lead_real <- lead_rows("genotype", genes_real) %>%
  select(chr, gene_name, test, variant_id, all_of(paste0("b_", CELLTYPES))) %>%
  inner_join(bulk %>% filter(bulk_qval <= 0.05), by = "gene_name")
lead_real$r <- NA_real_
for (ch in unique(lead_real$chr)) {
  pre <- file.path(out_dir, "genotypes", ch)
  bim <- fread(paste0(pre, ".bim"), header = FALSE)
  rows <- which(lead_real$chr == ch)
  k <- which(bim$V2 %in% unique(c(lead_real$variant_id[rows], lead_real$bulk_variant[rows])))
  pl <- read.plink(paste0(pre, ".bed"), paste0(pre, ".bim"), paste0(pre, ".fam"), select.snps = k)
  G <- as(pl$genotypes, "numeric")
  alt <- sapply(strsplit(colnames(G), "_"), `[`, 4)
  flip <- pl$map$allele.2 != alt
  G[, flip] <- 2 - G[, flip]
  for (i in rows) if (all(c(lead_real$variant_id[i], lead_real$bulk_variant[i]) %in% colnames(G)))
    lead_real$r[i] <- cor(G[, lead_real$variant_id[i]], G[, lead_real$bulk_variant[i]], use = "complete.obs")
}
same_signal <- lead_real %>% filter(r^2 >= 0.8) %>% mutate(bulk_slope_oriented = bulk_slope * sign(r))
bulk_tab <- bind_rows(lapply(CELLTYPES, function(ct) {
  d <- same_signal %>% filter(test == ct)
  data.frame(celltype = ct, egenes_same_signal_as_bulk = nrow(d),
             rho_effect_vs_bulk = if (nrow(d) >= 10) cor(d[[paste0("b_", ct)]], d$bulk_slope_oriented, method = "spearman") else NA,
             same_sign_as_bulk = mean(sign(d[[paste0("b_", ct)]]) == sign(d$bulk_slope_oriented)))
}))
cat("\nB. each cell type's effect vs the bulk eQTL slope (same signal, r2 >= 0.8):\n")
print(bulk_tab, digits = 2)
pB <- ggplot(same_signal %>% mutate(celltype = factor(test, levels = CELLTYPES),
                                    effect = case_when(test == "Neutrophils" ~ b_Neutrophils, test == "Monocytes" ~ b_Monocytes,
                                                       test == "Bcells" ~ b_Bcells, test == "CD4" ~ b_CD4,
                                                       test == "CD8" ~ b_CD8, test == "NK" ~ b_NK)),
             aes(bulk_slope_oriented, pmax(pmin(effect, 4), -4), colour = celltype)) +
  geom_hline(yintercept = 0, colour = "grey70") + geom_vline(xintercept = 0, colour = "grey70") +
  geom_point(size = 0.8, alpha = 0.6) +
  facet_wrap(~ celltype, nrow = 1) +
  scale_colour_manual(values = CT_COLS, guide = "none") +
  labs(x = "GTEx bulk whole-blood eQTL slope", y = "cell-type effect (clipped at +-4)",
       title = "B. Cell-type effect vs bulk effect (eGenes with the same lead signal as bulk)")

## ================= C. discovery overlap with bulk =================
tested <- unique(top$gene_name[top$set == "genotype"])
bulk_egenes <- intersect(bulk$gene_name[bulk$bulk_qval <= 0.05], tested)
# UpSet: the full bulk set (also the ~9,300 bulk-only eGenes); bar heights on a log10 scale so small groups stay visible
sets <- c(genes_real, list(`bulk (GTEx)` = bulk_egenes))
overlap_tab <- bind_rows(lapply(CELLTYPES, function(ct) {
  g <- genes_real[[ct]]
  data.frame(celltype = ct, egenes = length(g), also_bulk_egene = sum(g %in% bulk_egenes),
             share_bulk = mean(g %in% bulk_egenes),
             only_this_celltype = sum(!g %in% unlist(genes_real[setdiff(CELLTYPES, ct)])),
             new_vs_bulk = sum(!g %in% bulk_egenes))
})) %>% mutate(share_bulk_all_tested_genes = mean(tested %in% bulk_egenes))
cat("\nC. overlap with GTEx bulk whole-blood eGenes:\n")
print(overlap_tab, digits = 2)
upset_file <- file.path(fig_dir, "07d_upset_bulk.png")
# one row per gene, TRUE/FALSE per set; every bar labelled with its count (so the small groups stay readable next to
# the ~9,300 bulk-only genes)
membership <- data.frame(gene = unique(unlist(sets)))
for (s in names(sets)) membership[[s]] <- membership$gene %in% sets[[s]]
p_up <- ComplexUpset::upset(
  membership, intersect = rev(names(sets)), name = NULL, n_intersections = 30, sort_sets = FALSE,
  base_annotations = list(`genes in the intersection` =
    ComplexUpset::intersection_size(text = list(size = 2.6))),
  set_sizes = ComplexUpset::upset_set_size() + ylab("eGenes")) +
  labs(title = "Our cell-type eGenes and GTEx bulk whole-blood eGenes (all shown)")
ggsave(upset_file, p_up, width = 14, height = 7, dpi = 200)
cat("figure:", upset_file, "\n")

## ================= D. the OneK1K reference: is the sharing structured like real biology? =================
# leakage = every pair alike and high; noise = every pair alike and ~0; biology (OneK1K, single cells): related cell
# types share more (CD4-CD8 > T-NK > lymphoid-monocyte). OneK1K matrix, built like ours: genes that are eGenes
# (q < 0.05) in either type and have the same lead variant in both -> Spearman rho of the two slopes.
# (no neutrophils in OneK1K, a PBMC study -> 5 cell types, 10 pairs)
onek1k_dir <- "/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/rerun_LM22_aaron_1103/external/onek1k_top_eqtl/unzipped/OneK1K_TensorQTL_top_eQTL_summary"
OK_MAP <- c(Monocytes = "Mono_C", Bcells = "B_IN", CD4 = "CD4_NC", CD8 = "CD8_NC", NK = "NK")
ok <- lapply(OK_MAP, function(ty) bind_rows(lapply(CHRS, function(ch)
  fread(file.path(onek1k_dir, paste0("OneK1K_", ty, ".sig_cis_qtl_pairs.", ch, ".csv")),
        select = c("phenotype_id", "variant_id", "slope", "qval")))))
CT5 <- names(OK_MAP)
ok_pairs <- bind_rows(lapply(seq_along(CT5), function(i) bind_rows(lapply(CT5[-(1:i)], function(b) {
  a <- CT5[i]
  m <- inner_join(ok[[a]], ok[[b]], by = c("phenotype_id", "variant_id"), suffix = c("_a", "_b")) %>%
    filter(qval_a < 0.05 | qval_b < 0.05)
  data.frame(celltype_a = a, celltype_b = b, onek1k_genes = nrow(m),
             onek1k_rho = cor(m$slope_a, m$slope_b, method = "spearman"),
             onek1k_share_a_egenes_also_b = mean(ok[[b]]$qval[match(ok[[a]]$phenotype_id[ok[[a]]$qval < 0.05], ok[[b]]$phenotype_id)] < 0.05, na.rm = TRUE))
}))))
# ours, the same 10 pairs: rho above the shuffled baseline (the joint estimation alone gives a negative rho)
ours10 <- pairs %>%
  filter(as.character(celltype_a) %in% CT5, as.character(celltype_b) %in% CT5) %>%
  mutate(celltype_a = as.character(celltype_a), celltype_b = as.character(celltype_b),
         key = paste(pmin(celltype_a, celltype_b), pmax(celltype_a, celltype_b))) %>%
  group_by(key, data) %>% summarise(rho = mean(rho, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = data, values_from = rho) %>%
  mutate(our_rho_above_shuffled = `real genotypes` - `shuffled genotypes (mean of 5)`)
compare <- ok_pairs %>%
  mutate(key = paste(pmin(celltype_a, celltype_b), pmax(celltype_a, celltype_b))) %>%
  left_join(ours10, by = "key") %>%
  rename(our_rho = `real genotypes`, our_rho_shuffled = `shuffled genotypes (mean of 5)`) %>%
  mutate(pair = paste(celltype_a, "-", celltype_b)) %>%
  arrange(desc(onek1k_rho))
rank_agreement <- cor.test(compare$onek1k_rho, compare$our_rho_above_shuffled, method = "spearman", exact = FALSE)
cat("\nD. pair-by-pair: OneK1K effect correlation vs ours (above the shuffled baseline):\n")
print(compare %>% select(pair, onek1k_genes, onek1k_rho, our_rho, our_rho_shuffled, our_rho_above_shuffled) %>% as.data.frame(), digits = 2)
cat("rank agreement of the 10 pairs (Spearman):", round(rank_agreement$estimate, 2), " p =", signif(rank_agreement$p.value, 2), "\n")
pD <- ggplot(compare, aes(onek1k_rho, our_rho_above_shuffled, label = pair)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey60", linetype = "dashed") +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(size = 3.5) +
  labs(x = "OneK1K (single cells): effect correlation between the two cell types",
       y = "ours: effect correlation above the shuffled baseline",
       title = "D. Is the sharing structured like real biology? (each point = a pair of cell types)",
       subtitle = sprintf("rank agreement of the 10 pairs: Spearman %.2f (p = %.2g). Leakage would put every pair at the same high value.",
                          rank_agreement$estimate, rank_agreement$p.value))
d_file <- file.path(fig_dir, "07d_onek1k_structure.png")
ggsave(d_file, pD, width = 9, height = 7, dpi = 200)
cat("figure:", d_file, "\n")
fwrite(compare, file.path(out_dir, paste0("07d_leakage_onek1k_structure", SCAN_TAG, ".csv")))

## ================= save =================
fig_file <- file.path(fig_dir, "07d_leakage.png")
ggsave(fig_file, (pA | pA2) / pB + plot_layout(heights = c(1.3, 1)), width = 18, height = 12, dpi = 200)
cat("figure:", fig_file, "\n")
fwrite(pairs, file.path(out_dir, paste0("07d_leakage_pairs", SCAN_TAG, ".csv")))
fwrite(bulk_tab, file.path(out_dir, paste0("07d_leakage_bulk", SCAN_TAG, ".csv")))
fwrite(overlap_tab, file.path(out_dir, paste0("07d_leakage_overlap", SCAN_TAG, ".csv")))
message("07d_leakage_check done")
