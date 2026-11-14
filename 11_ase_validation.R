#author: haim krupkin+claude
#date: 11.12.2026
#description: step 11. Allele-specific expression (ASE) validation of the cell-type eGenes, as in Kim-Hellmuth et al.
# 2020 (Science, Fig. 2B). This check does NOT use the deconvolution model: within one person, the two copies of a gene
# are in the same cells, so the ratio between the ALT and REF copy (allelic fold change, aFC) is free of composition
# and batch confounding between people. If an eQTL is active in cell type c, heterozygotes with more of c should show
# a stronger allelic imbalance: aFC correlates with the proportion of c, in the direction of our effect.
# Per eGene (calibrated hit, step 9 hit list) at its lead variant:
#   1. heterozygotes at the lead variant (population-phased WGS genotype, phASER PW field)
#   2. their haplotype read counts for the gene in whole blood (phASER, phase kept across genes = same haplotypes
#      as PW); aFC = log2((ALT-haplotype reads + 0.5) / (REF-haplotype reads + 0.5)); >= 8 reads per sample
#   3. drop outlying aFC (Hampel: |aFC - median| >= 4.5 x MAD), as the paper; need >= 10 heterozygotes
#   4. Pearson correlation of aFC with the cell type's proportion; validated = p < 0.05 and the sign of our effect
#   Baseline: the same test with the proportions shuffled between people (what chance gives).
# Inputs : figures/09_3_all_hits.csv, 06_scan_<chr>_top_variants.csv, 05_proportions_6classes.csv, 02_covariates.csv,
#          GTEx v8 phASER (gw_phased matrix + merged VCF), GENCODE v26
# Outputs: output/11_ase_validation_genes.csv (per eGene), output/11_ase_validation_summary.csv (per cell type)
suppressMessages(library(dplyr))
suppressMessages(library(data.table))

## ================= settings =================
pipe_dir <- "/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/for_evin_cibersorting_expression/rerun_LM22_aaron_1103/final_pipeline"
out_dir  <- file.path(pipe_dir, "output")
tmp_dir  <- file.path(out_dir, "11_tmp")
dir.create(tmp_dir, showWarnings = FALSE)
phaser   <- "/labs/smontgom/grps_smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/ase/GTEx_Analysis_v8_phASER"
gtf_file <- "/labs/smontgom/grps_smontgom/shared/GTEx/all_data/GTEx_Analysis_2017-06-05_v8/references/gencode.v26.GRCh38.genes.gtf"
TABIX    <- "/home/hkrupkin/envs/catca-edit/bin/tabix"
CELLTYPES <- c("Neutrophils", "Monocytes", "Bcells", "CD4", "CD8", "NK")
MIN_READS <- 8
MIN_HETS  <- 10
SCAN_TAG  <- Sys.getenv("SCAN_TAG", "_epc20")  # step-6 OUT_TAG (hit list from output/figures<SCAN_TAG>/). 11.13: main = 20 hidden-noise PCs; "" = without
set.seed(1)

## ================= the eGenes and their lead variants =================
hits <- read.csv(file.path(out_dir, paste0("figures", SCAN_TAG), "09_3_all_hits.csv"))
top <- bind_rows(lapply(unique(hits$chr), function(ch) {
  fread(file.path(out_dir, paste0("06_scan_", ch, "_top_variants", SCAN_TAG, ".csv"))) %>%
    filter(set == "genotype") %>%
    select(gene_name, test, variant_id, all_of(paste0("b_", CELLTYPES)))
}))
hits <- hits %>%
  left_join(top, by = c("gene_name", "celltype" = "test")) %>%
  mutate(our_b = case_when(celltype == "Neutrophils" ~ b_Neutrophils, celltype == "Monocytes" ~ b_Monocytes,
                           celltype == "Bcells" ~ b_Bcells, celltype == "CD4" ~ b_CD4,
                           celltype == "CD8" ~ b_CD8, celltype == "NK" ~ b_NK)) %>%
  select(chr, gene_name, celltype, variant_id, our_b)
stopifnot(!anyNA(hits$variant_id))

# gene symbol -> GENCODE v26 gene ID (the phASER matrix uses the same gene models)
gtf <- fread(cmd = paste("grep -P '\\tgene\\t'", gtf_file), header = FALSE)
gtf$gene_name <- sub('.*gene_name "([^"]+)".*', "\\1", gtf$V9)
gtf$gene_id <- sub('.*gene_id "([^"]+)".*', "\\1", gtf$V9)
hits$gene_id <- gtf$gene_id[match(hits$gene_name, gtf$gene_name)]

## ================= whole-blood samples and proportions =================
cov_df <- read.csv(file.path(out_dir, "02_covariates.csv"))
W_df <- read.csv(file.path(out_dir, "05_proportions_6classes.csv"))
W_df$SUBJID <- cov_df$SUBJID[match(W_df$SAMPID, cov_df$SAMPID)]

## ================= haplotype counts of the eGenes in whole blood (phASER matrix) =================
# streamed with awk: keep the header + the eGene rows, and only the whole-blood columns
counts_file <- file.path(tmp_dir, "phaser_counts_hits.tsv")
writeLines(unique(hits$gene_id), file.path(tmp_dir, "gene_ids.txt"))
writeLines(W_df$SAMPID, file.path(tmp_dir, "sample_ids.txt"))
awk_prog <- paste0(
  "FILENAME==ARGV[1]{g[$1]=1; next} FILENAME==ARGV[2]{s[$1]=1; next} ",
  "FNR==1{n=0; for(i=1;i<=NF;i++) if(i==2 || ($i in s)) keep[++n]=i} ",
  "FNR==1 || ($2 in g){line=$(keep[1]); for(j=2;j<=n;j++) line=line\"\\t\"$(keep[j]); print line}")
cmd <- paste0("zcat ", file.path(phaser, "phASER_GTEx_v8_matrix.gw_phased.txt.gz"),
              " | awk -F'\\t' '", awk_prog, "' ", file.path(tmp_dir, "gene_ids.txt"), " ",
              file.path(tmp_dir, "sample_ids.txt"), " - > ", counts_file)
system(cmd)
counts <- fread(counts_file)
cat("eGenes with phASER counts:", nrow(counts), "of", length(unique(hits$gene_id)),
    "; whole-blood samples:", ncol(counts) - 1, "\n")

## ================= phased genotypes at the lead variants (phASER VCF, PW = population phase) =================
regions <- hits %>% distinct(variant_id) %>%
  mutate(chr = sub("_.*", "", variant_id), pos = as.integer(sapply(strsplit(variant_id, "_"), `[`, 2))) %>%
  arrange(chr, pos)
fwrite(regions %>% select(chr, pos), file.path(tmp_dir, "regions.tsv"), sep = "\t", col.names = FALSE)
vcf <- file.path(phaser, "phASER_GTEx_v8_merged.vcf.gz")
hdr <- system(paste(TABIX, "-H", vcf, "| tail -1"), intern = TRUE)
vcf_samples <- strsplit(hdr, "\t")[[1]][-(1:9)]
system(paste(TABIX, "-R", file.path(tmp_dir, "regions.tsv"), vcf, ">", file.path(tmp_dir, "lead_variants.vcf")))
geno <- fread(file.path(tmp_dir, "lead_variants.vcf"), header = FALSE, sep = "\t")
setnames(geno, c("chrom", "pos", "id", "ref", "alt", "qual", "filter", "info", "format", vcf_samples))
geno <- geno[id %in% regions$variant_id]
pw_index <- which(strsplit(geno$format[1], ":")[[1]] == "PW")
stopifnot(all(geno$format == geno$format[1]))
cat("lead variants found in the phASER VCF:", nrow(geno), "of", nrow(regions), "\n")

## ================= per eGene: aFC of heterozygotes vs the cell type's proportion =================
# blood sample of each subject (one per subject in our 670)
subj <- W_df$SUBJID
test_one <- function(k) {
  h <- hits[k, ]
  gi <- which(geno$id == h$variant_id)
  ci <- which(counts$name == h$gene_id)
  if (length(gi) != 1 || length(ci) != 1) return(NULL)
  pw <- sapply(strsplit(unlist(geno[gi, subj, with = FALSE]), ":"), `[`, pw_index)   # e.g. "0|1"
  hap <- strsplit(unlist(counts[ci, W_df$SAMPID, with = FALSE]), "|", fixed = TRUE)
  a <- as.numeric(sapply(hap, `[`, 1))
  b <- as.numeric(sapply(hap, `[`, 2))
  # ALT haplotype reads / REF haplotype reads: "0|1" -> ALT on haplotype B; "1|0" -> ALT on haplotype A
  alt_reads <- ifelse(pw == "0|1", b, ifelse(pw == "1|0", a, NA))
  ref_reads <- ifelse(pw == "0|1", a, ifelse(pw == "1|0", b, NA))
  d <- data.frame(alt = alt_reads, ref = ref_reads, W_df[, CELLTYPES]) %>%
    filter(!is.na(alt), alt + ref >= MIN_READS) %>%
    mutate(afc = log2((alt + 0.5) / (ref + 0.5)))
  d <- d %>% filter(abs(afc - median(afc)) < 4.5 * mad(afc, constant = 1))     # Hampel outliers out
  if (nrow(d) < MIN_HETS || sd(d$afc) == 0) return(NULL)
  x <- d[[h$celltype]]
  ct <- cor.test(d$afc, x)
  x_shuffled <- sample(x)                                                        # baseline: chance
  cs <- cor.test(d$afc, x_shuffled)
  # also: the correlation with every cell type (is it strongest for the called one?)
  r_all <- sapply(CELLTYPES, function(c2) cor(d$afc, d[[c2]]))
  # is there allelic imbalance at all? (one-sample t-test of aFC; the paper kept only genes with ASE aFC P < 0.05)
  p_imbalance <- t.test(d$afc)$p.value
  data.frame(h, n_hets = nrow(d), median_afc = median(d$afc), p_imbalance = p_imbalance,
             r = unname(ct$estimate), p = ct$p.value,
             r_shuffled = unname(cs$estimate), p_shuffled = cs$p.value,
             strongest_celltype = CELLTYPES[which.max(abs(r_all))],
             t(setNames(r_all, paste0("r_", CELLTYPES))))
}
res <- bind_rows(lapply(seq_len(nrow(hits)), test_one))
res <- res %>%
  mutate(validated = p < 0.05 & sign(r) == sign(our_b),
         validated_shuffled = p_shuffled < 0.05 & sign(r_shuffled) == sign(our_b))

## ================= summary per cell type =================
summary_tab <- res %>%
  group_by(celltype) %>%
  # (shares first: in summarise a new column with the same name would replace "validated")
  summarise(egenes_testable = n(),
            median_hets = median(n_hets),
            share_validated = mean(validated),
            share_validated_shuffled = mean(validated_shuffled),
            n_validated = sum(validated),
            share_p05_any_sign = mean(p < 0.05),
            # as the paper: only eGenes with allelic imbalance (ASE aFC P < 0.05)
            with_imbalance = sum(p_imbalance < 0.05),
            share_validated_with_imbalance = mean(validated[p_imbalance < 0.05]),
            share_validated_shuffled_with_imbalance = mean(validated_shuffled[p_imbalance < 0.05]),
            share_strongest_is_called = mean(strongest_celltype == celltype),
            .groups = "drop") %>%
  right_join(hits %>% count(celltype, name = "egenes"), by = "celltype") %>%
  relocate(egenes, .after = celltype)
summary_tab$celltype <- factor(summary_tab$celltype, levels = CELLTYPES)
summary_tab <- summary_tab %>% arrange(celltype)
print(as.data.frame(summary_tab), digits = 2)

fwrite(res, file.path(out_dir, paste0("11_ase_validation_genes", SCAN_TAG, ".csv")))
fwrite(summary_tab, file.path(out_dir, paste0("11_ase_validation_summary", SCAN_TAG, ".csv")))
cat("outputs:", file.path(out_dir, paste0("11_ase_validation_genes", SCAN_TAG, ".csv")), file.path(out_dir, paste0("11_ase_validation_summary", SCAN_TAG, ".csv")), "\n")
message("11_ase_validation done")
