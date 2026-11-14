#author: haim krupkin+claude
#date: 11.14.2026
#description: step 12a. Candidate variants for the cell-type trans-eQTL scan (step 12).
# A genome-wide trans scan (every variant x every gene) is ~1,000x the cis scan and has little power at n = 670, so
# (as eQTLGen and most trans studies) only candidate variants are tested against all genes:
#   CAND = cis_leads (default): the lead variant of every calibrated cell-type cis-eGene of step 6/7 (empirical FDR < 5%,
#          same calibration as step 7). Trans effects are often carried by a cis effect on a regulator, so these are
#          the variants most likely to have one.
#   CAND = <path to a file with a column variant_id> (e.g. GWAS blood-trait lead variants in GTEx id format
#          chr1_12345_A_G_b38): those variants instead (source = "file").
# Variants in strong LD would test the same signal twice: within 1 Mb, a variant with r2 >= 0.8 (GTEx genotypes) to a
# variant already kept is dropped (the one from more cis-eGenes / smaller p is kept first).
# Inputs : output/06_scan_<chr>_top_variants<SCAN_TAG>.csv (all CHRS), output/genotypes/<chr>.{bed,bim,fam}
# Output : output/12a_trans_candidates<SCAN_TAG>.csv  (variant_id, chr, pos, source cell types, cis genes, best cis p)
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(snpStats))

## ================= settings =================
pipe_dir <- "/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/rerun_LM22_aaron_1103/final_pipeline"
out_dir  <- file.path(pipe_dir, "output")
CHRS     <- strsplit(Sys.getenv("CHRS", paste0("chr", 1:22, collapse = " ")), " ")[[1]]
SCAN_TAG <- Sys.getenv("SCAN_TAG", "_epc20")
CAND     <- Sys.getenv("CAND", "cis_leads")
CELLTYPES <- c("Neutrophils", "Monocytes", "Bcells", "CD4", "CD8", "NK")
FDR      <- 0.05
LD_R2    <- 0.8      # pruning threshold
LD_WIN   <- 1e6      # pruning window (bp)

## ================= candidates =================
if (CAND == "cis_leads") {
  # main scans of all chromosomes (real + shuffled), the 6 cell-type tests only
  main <- bind_rows(lapply(CHRS, function(ch) {
    fread(file.path(out_dir, paste0("06_scan_", ch, "_top_variants", SCAN_TAG, ".csv")),
          select = c("gene_name", "variant_id", "test", "set", "p_min", "p_gene", "chr"))
  }))
  main <- main %>% filter(test %in% CELLTYPES)
  # empirical p against the pooled shuffles, then BH (the same calibration as step 7)
  empirical_fdr <- function(p_obs, p_null) {
    p_null <- sort(p_null)
    n_le <- findInterval(p_obs, p_null)                       # number of null values <= p_obs
    p.adjust((1 + n_le) / (1 + length(p_null)), "BH")
  }
  real <- main %>% filter(set == "genotype")
  real$fdr_emp <- NA_real_
  for (tt in CELLTYPES) {
    null_p <- main$p_gene[main$test == tt & main$set != "genotype"]
    is_tt <- real$test == tt
    real$fdr_emp[is_tt] <- empirical_fdr(real$p_gene[is_tt], null_p)
  }
  hits <- real %>% filter(fdr_emp < FDR)
  print(table(hits$test))
  # one row per variant: which cell types / genes it is a cis lead for
  cand <- hits %>%
    group_by(variant_id, chr) %>%
    summarise(source = paste(sort(unique(test)), collapse = ";"),
              cis_genes = paste(sort(unique(gene_name)), collapse = ";"),
              n_cis = n_distinct(gene_name, test),
              best_cis_p = min(p_min),
              .groups = "drop")
} else {
  # a user-supplied list of variants
  cand <- fread(CAND) %>%
    distinct(variant_id) %>%
    mutate(chr = sub("_.*", "", variant_id), source = "file", cis_genes = NA_character_, n_cis = 0L, best_cis_p = NA_real_)
}
cand$pos <- as.numeric(sapply(strsplit(cand$variant_id, "_"), `[`, 2))
cat("candidate variants before LD pruning:", nrow(cand), "\n")

## ================= LD pruning, per chromosome =================
# order: variants leading more cis-eGenes first, then by smaller cis p
cand <- cand %>% arrange(desc(n_cis), best_cis_p)
prune_chr <- function(ch) {
  d <- cand %>% filter(chr == ch)
  if (nrow(d) <= 1) return(d)
  pfx <- file.path(out_dir, "genotypes", ch)
  # keep the candidates present in the genotype file (variant ids = column 2 of the bim)
  bim_ids <- fread(paste0(pfx, ".bim"), header = FALSE, select = 2)$V2
  found <- d$variant_id %in% bim_ids
  if (any(!found)) cat(ch, ": ", sum(!found), " candidate variants not in the genotype file (dropped)\n", sep = "")
  d <- d[found, ]
  if (nrow(d) <= 1) return(d)
  pl <- read.plink(paste0(pfx, ".bed"), paste0(pfx, ".bim"), paste0(pfx, ".fam"), select.snps = d$variant_id)
  G <- as(pl$genotypes[, d$variant_id], "numeric")
  kept <- logical(nrow(d))
  for (i in seq_len(nrow(d))) {
    near_kept <- which(kept & abs(d$pos - d$pos[i]) <= LD_WIN)
    if (length(near_kept) > 0) {
      r2 <- cor(G[, i], G[, near_kept, drop = FALSE], use = "pairwise.complete.obs")^2
      if (any(r2 >= LD_R2, na.rm = TRUE)) next
    }
    kept[i] <- TRUE
  }
  d[kept, ]
}
cand <- bind_rows(lapply(unique(cand$chr), prune_chr))
cat("candidate variants after LD pruning (r2 <", LD_R2, "within", LD_WIN / 1e6, "Mb):", nrow(cand), "\n")
print(table(cand$chr))

## ================= save =================
cand_tag <- if (CAND == "cis_leads") "" else paste0("_", tools::file_path_sans_ext(basename(CAND)))
out_file <- file.path(out_dir, paste0("12a_trans_candidates", cand_tag, SCAN_TAG, ".csv"))
write.csv(cand %>% select(variant_id, chr, pos, source, cis_genes, n_cis, best_cis_p), out_file, row.names = FALSE)
out_file
message("12a_trans_candidates done")
