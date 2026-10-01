#!/usr/bin/env Rscript
# 02_host_filter.R
# Build a host-free bacterial SQM object:
#   1. flag contigs whose sequence aligns to the Hydra 105 genome (PAF from 01_contigs_vs_hydra.sh)
#   2. keep contigs classified as Bacteria that are NOT host-flagged (TPM rescaled within the subset)
#   3. REPORT (do not drop) Eukaryota-labelled ORFs on the kept bacterial contigs
#   4. REPORT (do not drop) kept bacterial contigs that sit in host bins (metabat2.4, metabat2.2)
#
# Input (read-only):
#   /mnt/hdd3/sqm_data_cache.rds                         full loadSQM() object (~4.8 GB RAM)
#   meta/result/host_filter/contigs_vs_hydra105.paf
# Output (meta/result/host_filter/):
#   sqm_bacteria_hostfree.rds        SQM object used for all downstream bacterial analyses
#   contig_host_flags.tsv.gz         per-contig taxonomy, GC, length, bin, Hydra coverage, flags
#   reads_by_category.tsv            per-sample reads in each contig category (long format)
#   sample_summary.tsv               per-sample bacterial reads before/after + Eukaryota-ORF share
#   eukaryota_orfs_on_bacterial_contigs.tsv
#   bacterial_contigs_in_host_bins.tsv
#   eukaryote_marker_pfam_check.tsv
#   host_filter_summary.txt
#
# Usage: Rscript meta/scripts/02_host_filter.R

suppressPackageStartupMessages({
  library(SQMtools)
  library(data.table)
})

cache_rds <- "/mnt/hdd3/sqm_data_cache.rds"
out_dir   <- path.expand("~/HVBMT_project/meta/result/host_filter")
paf_file  <- file.path(out_dir, "contigs_vs_hydra105.paf")

# Host-flag thresholds (WORKFLOW.md, 2026-09-30). Cultures are strain 105, same as the reference.
min_identity <- 0.90   # per alignment: matching bases / alignment block length
min_coverage <- 0.50   # per contig: fraction of length covered by passing alignments
host_bins    <- c("sqm_noninduced.metabat2.4", "sqm_noninduced.metabat2.2")

# Domains expected only in eukaryotes; should be ~0 in the host-free bacterial subset
marker_pfams <- c(PF00001 = "7tm_1 GPCR (rhodopsin family)",
                  PF00002 = "7tm_2 GPCR (secretin family)",
                  PF00003 = "7tm_3 GPCR (class C)",
                  PF00125 = "Core histone",
                  PF00022 = "Actin")

stopifnot(file.exists(cache_rds), file.exists(paf_file), file.size(paf_file) > 0)

# SQM count matrices can carry "Raw read count <sample>" column names; normalise to sample IDs
clean_cols <- function(m) {
  colnames(m) <- sub("^Raw read count ", "", colnames(m))
  m
}

# ---- 1. Hydra coverage per contig from the PAF ---------------------------------------------
message("Reading PAF ...")
paf <- fread(cmd = sprintf("cut -f1-11 '%s'", paf_file), sep = "\t", header = FALSE,
             col.names = c("contig", "qlen", "qstart", "qend", "strand", "target",
                           "tlen", "tstart", "tend", "nmatch", "alnlen"))
paf[, identity := nmatch / alnlen]

# Union of passing alignment intervals on each contig, so overlapping hits are not double-counted
hits <- paf[identity >= min_identity]
setorder(hits, contig, qstart)
hits[, prev_end := shift(cummax(qend), fill = -1L), by = contig]
hits[, block := cumsum(qstart > prev_end), by = contig]
host_cov <- hits[, .(s = min(qstart), e = max(qend)), by = .(contig, qlen, block)][
  , .(hydra_cov = sum(e - s) / qlen[1]), by = contig]
best_id <- paf[, .(hydra_best_identity = max(identity)), by = contig]

# ---- 2. Contig table -------------------------------------------------------------------------
message("Loading SQM cache (~4.8 GB) ...")
SQM <- readRDS(cache_rds)
stopifnot(inherits(SQM, "SQM"))

contig_ids <- rownames(SQM$contigs$table)
missing_in_sqm <- setdiff(unique(paf$contig), contig_ids)
if (length(missing_in_sqm)) {
  stop(sprintf("%d PAF contig IDs are not in the SQM object (e.g. %s)",
               length(missing_in_sqm), missing_in_sqm[1]))
}

ctg <- data.table(
  contig       = contig_ids,
  superkingdom = SQM$contigs$tax[contig_ids, "superkingdom"],
  length       = SQM$contigs$table[contig_ids, "Length"],
  gc           = SQM$contigs$table[contig_ids, "GC perc"],
  bin          = if (is.null(SQM$contigs$bins)) "No bin" else SQM$contigs$bins[contig_ids, 1]
)
ctg <- merge(ctg, host_cov, by = "contig", all.x = TRUE, sort = FALSE)
ctg <- merge(ctg, best_id,  by = "contig", all.x = TRUE, sort = FALSE)
ctg[is.na(hydra_cov), hydra_cov := 0]
ctg[, host_aligned := hydra_cov >= min_coverage]
ctg[, in_host_bin  := bin %in% host_bins]
ctg[, category := fcase(
  superkingdom == "Bacteria" & !host_aligned, "Bacteria (kept)",
  superkingdom == "Bacteria" &  host_aligned, "Bacteria (host-aligned, removed)",
  superkingdom == "Eukaryota",                "Eukaryota",
  superkingdom == "Unclassified" & host_aligned, "Unclassified (host-aligned)",
  superkingdom == "Unclassified",             "Unclassified (no host hit)",
  default = "Other (Archaea, Viruses, No CDS)")]

# ---- 3. Reads per category per sample --------------------------------------------------------
contig_abund <- clean_cols(SQM$contigs$abund)[ctg$contig, , drop = FALSE]
total_reads  <- SQM$total_reads
stopifnot(setequal(colnames(contig_abund), names(total_reads)))
total_reads  <- total_reads[colnames(contig_abund)]

cat_reads <- rowsum(contig_abund, ctg$category)
cat_reads <- rbind(cat_reads, Unmapped = total_reads - colSums(contig_abund))
reads_long <- as.data.table(as.table(cat_reads))
setnames(reads_long, c("category", "sample", "reads"))
reads_long[, pct := 100 * reads / total_reads[as.character(sample)]]

# ---- 4. Host-free bacterial subset -----------------------------------------------------------
keep_contigs <- ctg[category == "Bacteria (kept)", contig]
message(sprintf("Subsetting to %d bacterial contigs ...", length(keep_contigs)))
SQM_bac <- subsetContigs(SQM, keep_contigs, rescale_tpm = TRUE, rescale_copy_number = TRUE,
                         recalculate_bin_stats = FALSE)

# ---- 5. Report Eukaryota-labelled ORFs on kept bacterial contigs (not removed) --------------
orf_abund <- clean_cols(SQM_bac$orfs$abund)
orf_sk    <- SQM_bac$orfs$tax[, "superkingdom"]
euk_orfs  <- rownames(SQM_bac$orfs$tax)[orf_sk == "Eukaryota"]

ann_cols <- intersect(c("Contig ID", "Length NT", "GC perc", "Gene name", "KEGG ID", "KEGGFUN",
                        "COG ID", "COGFUN", "PFAM"), colnames(SQM_bac$orfs$table))
euk_tab <- as.data.table(SQM_bac$orfs$table[euk_orfs, ann_cols, drop = FALSE],
                         keep.rownames = "orf")
if (length(euk_orfs)) {
  lineage <- function(tax) apply(tax, 1, paste, collapse = ";")
  euk_tab[, orf_tax    := lineage(SQM_bac$orfs$tax[euk_orfs, , drop = FALSE])]
  euk_tab[, contig_tax := lineage(SQM_bac$contigs$tax[`Contig ID`, , drop = FALSE])]
  euk_tab[, total_reads := rowSums(orf_abund[euk_orfs, , drop = FALSE])]
  euk_tab <- cbind(euk_tab, orf_abund[euk_orfs, , drop = FALSE])
  setorder(euk_tab, -total_reads)
}

# ---- 6. Report kept bacterial contigs that sit in host bins (not removed) -------------------
host_bin_tab <- ctg[category == "Bacteria (kept)" & in_host_bin,
                   .(contig, bin, length, gc, hydra_cov, hydra_best_identity)]
host_bin_tab[, total_reads := rowSums(contig_abund[contig, , drop = FALSE])]
setorder(host_bin_tab, -total_reads)

# ---- 7. Eukaryote marker PFAM check ----------------------------------------------------------
pfam_reads <- function(m, id) {
  r <- grep(paste0("^", id), rownames(m))
  if (length(r)) sum(m[r, ]) else 0
}
pfam_check <- data.table(
  pfam             = names(marker_pfams),
  description      = unname(marker_pfams),
  reads_all_contigs = sapply(names(marker_pfams), pfam_reads, m = SQM$functions$PFAM$abund),
  reads_bacteria_hostfree = sapply(names(marker_pfams), pfam_reads,
                                   m = SQM_bac$functions$PFAM$abund))

# ---- 8. Per-sample summary -------------------------------------------------------------------
get_reads <- function(cat) {
  if (cat %in% rownames(cat_reads)) cat_reads[cat, ] else setNames(numeric(ncol(cat_reads)), colnames(cat_reads))
}
bac_orf_reads <- colSums(orf_abund)
euk_orf_reads <- colSums(orf_abund[euk_orfs, , drop = FALSE])
sample_summary <- data.table(
  sample              = names(total_reads),
  total_reads         = total_reads,
  bacteria_before     = get_reads("Bacteria (kept)") + get_reads("Bacteria (host-aligned, removed)"),
  bacteria_removed    = get_reads("Bacteria (host-aligned, removed)"),
  bacteria_kept       = get_reads("Bacteria (kept)"))
sample_summary[, pct_bacteria_kept := 100 * bacteria_kept / total_reads]
sample_summary[, bacterial_orf_reads := bac_orf_reads[sample]]
sample_summary[, eukaryota_orf_reads := euk_orf_reads[sample]]
sample_summary[, pct_eukaryota_orf  := 100 * eukaryota_orf_reads / bacterial_orf_reads]

# ---- 9. Write outputs ------------------------------------------------------------------------
message("Writing outputs ...")
saveRDS(SQM_bac, file.path(out_dir, "sqm_bacteria_hostfree.rds"))
fwrite(ctg,            file.path(out_dir, "contig_host_flags.tsv.gz"), sep = "\t")
fwrite(reads_long,     file.path(out_dir, "reads_by_category.tsv"), sep = "\t")
fwrite(sample_summary, file.path(out_dir, "sample_summary.tsv"), sep = "\t")
fwrite(euk_tab,        file.path(out_dir, "eukaryota_orfs_on_bacterial_contigs.tsv"), sep = "\t")
fwrite(host_bin_tab,   file.path(out_dir, "bacterial_contigs_in_host_bins.tsv"), sep = "\t")
fwrite(pfam_check,     file.path(out_dir, "eukaryote_marker_pfam_check.tsv"), sep = "\t")

bac <- ctg[superkingdom == "Bacteria"]
unc <- ctg[superkingdom == "Unclassified"]
summary_lines <- c(
  sprintf("Host filter run: %s", format(Sys.time(), "%Y-%m-%d %H:%M")),
  sprintf("Thresholds: identity >= %.2f, contig coverage >= %.2f", min_identity, min_coverage),
  sprintf("Contigs with any Hydra alignment: %d of %d", nrow(best_id), nrow(ctg)),
  sprintf("Bacteria contigs: %d; host-aligned and removed: %d (%.1f%%); kept: %d",
          nrow(bac), sum(bac$host_aligned), 100 * mean(bac$host_aligned), sum(!bac$host_aligned)),
  sprintf("Unclassified contigs host-aligned: %d of %d (%.1f%%)",
          sum(unc$host_aligned), nrow(unc), 100 * mean(unc$host_aligned)),
  sprintf("Bacterial reads removed as host: %.2f%% of bacterial reads (all samples)",
          100 * sum(sample_summary$bacteria_removed) / sum(sample_summary$bacteria_before)),
  sprintf("Kept bacterial reads per sample: %.2f-%.2f%% of total reads",
          min(sample_summary$pct_bacteria_kept), max(sample_summary$pct_bacteria_kept)),
  sprintf("Eukaryota-labelled ORFs on kept bacterial contigs (reported, not removed): %d; %.2f-%.2f%% of bacterial ORF reads",
          length(euk_orfs), min(sample_summary$pct_eukaryota_orf), max(sample_summary$pct_eukaryota_orf)),
  sprintf("Kept bacterial contigs in host bins (reported, not removed): %d", nrow(host_bin_tab)),
  "Eukaryote marker PFAM reads (all contigs -> host-free bacteria):",
  sprintf("  %s %-32s %10.0f -> %.0f", pfam_check$pfam, pfam_check$description,
          pfam_check$reads_all_contigs, pfam_check$reads_bacteria_hostfree))
writeLines(summary_lines, file.path(out_dir, "host_filter_summary.txt"))
writeLines(summary_lines)
