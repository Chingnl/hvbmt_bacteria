#!/usr/bin/env Rscript
# 03_host_filter_qc.R
# QC figures for the host filter. Reads only the small tables written by 02_host_filter.R
# (not the 4.8 GB SQM cache).
#
# Input:  meta/result/host_filter/{reads_by_category,sample_summary,eukaryota_orfs_on_bacterial_contigs,
#                                  eukaryote_marker_pfam_check}.tsv, contig_host_flags.tsv.gz
#         meta/HVBMT_metadata_noninduced.csv
# Output: meta/result/figures/03_*.{pdf,png,svg}
#         meta/result/host_filter/eukaryota_orfs_top_functions.tsv
#
# Usage: Rscript meta/scripts/03_host_filter_qc.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

meta_dir <- path.expand("~/HVBMT_project/meta")
hf_dir   <- file.path(meta_dir, "result", "host_filter")
fig_dir  <- file.path(meta_dir, "result", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Shared style ----------------------------------------------------------------------------
# Categories of interest take the first validated categorical slots (blue, orange, aqua);
# host and background categories are neutral greys so the bacterial signal stands out.
col_kept    <- "#2a78d6"
col_removed <- "#eb6834"
col_unc_new <- "#1baf7a"
category_cols <- c(
  "Bacteria (kept)"                  = col_kept,
  "Bacteria (host-aligned, removed)" = col_removed,
  "Unclassified (no host hit)"       = col_unc_new,
  "Unclassified (host-aligned)"      = "#a8a7a2",
  "Eukaryota"                        = "#6f6e69",
  "Other (Archaea, Viruses, No CDS)" = "#c9c8c3",
  "Unmapped"                         = "#e4e3df")

theme_qc <- theme_minimal(base_size = 10) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_line(colour = "#e4e3df", size = 0.3),
        axis.text          = element_text(colour = "#52514e"),
        axis.title         = element_text(colour = "#52514e"),
        strip.text         = element_text(face = "bold", colour = "#0b0b0b"),
        legend.position    = "bottom",
        plot.title         = element_text(face = "bold", colour = "#0b0b0b"),
        plot.subtitle      = element_text(colour = "#52514e"))

save_fig <- function(p, name, width, height) {
  base <- file.path(fig_dir, name)
  ggsave(paste0(base, ".pdf"), p, width = width, height = height, device = cairo_pdf)
  ggsave(paste0(base, ".png"), p, width = width, height = height, dpi = 300, type = "cairo")
  ggsave(paste0(base, ".svg"), p, width = width, height = height, device = grDevices::svg)
  message("Saved ", base, ".{pdf,png,svg}")
}

# ---- Metadata: order samples by day, then temperature, then replicate -------------------------
meta <- fread(file.path(meta_dir, "HVBMT_metadata_noninduced.csv"))
setorder(meta, Time, Temp, Replicate)
meta[, day   := factor(paste0("Day ", Time), levels = paste0("Day ", sort(unique(Time))))]
meta[, label := paste0(Temp, "°C r", Replicate)]
meta[, SampleID := factor(SampleID, levels = SampleID)]
sample_labels <- setNames(meta$label, as.character(meta$SampleID))

add_meta <- function(dt) {
  dt <- merge(dt, meta[, .(sample = as.character(SampleID), day)], by = "sample")
  dt[, sample := factor(sample, levels = levels(meta$SampleID))]
  dt
}

# ---- Figure 1: read composition per sample ---------------------------------------------------
reads <- add_meta(fread(file.path(hf_dir, "reads_by_category.tsv")))
reads[, category := factor(category, levels = rev(names(category_cols)))]

p1 <- ggplot(reads, aes(sample, pct, fill = category)) +
  geom_col(width = 0.8, colour = "white", size = 0.2) +
  facet_grid(~ day, scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = sample_labels) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = category_cols, breaks = names(category_cols), name = NULL) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(title = "Read composition after host filtering",
       subtitle = "Contig categories from SqueezeMeta taxonomy + minimap2 alignment to H. vulgaris 105",
       x = NULL, y = "% of reads") +
  theme_qc + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
save_fig(p1, "03_host_filter_read_composition", width = 10, height = 5.5)

# ---- Figure 2: bacterial reads kept vs removed per sample ------------------------------------
smry <- add_meta(fread(file.path(hf_dir, "sample_summary.tsv")))
bac_long <- melt(smry, id.vars = c("sample", "day"),
                 measure.vars = c("bacteria_kept", "bacteria_removed"),
                 variable.name = "status", value.name = "reads")
bac_long[, status := factor(status, levels = c("bacteria_removed", "bacteria_kept"),
                            labels = c("Host-aligned, removed", "Kept"))]

p2 <- ggplot(bac_long, aes(sample, reads / 1e6, fill = status)) +
  geom_col(width = 0.8, colour = "white", size = 0.2) +
  facet_grid(~ day, scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = sample_labels) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c("Kept" = col_kept, "Host-aligned, removed" = col_removed),
                    breaks = c("Kept", "Host-aligned, removed"), name = NULL) +
  labs(title = "Bacterial reads per sample",
       subtitle = "Reads on Bacteria-classified contigs, split by Hydra alignment",
       x = NULL, y = "Reads (millions)") +
  theme_qc + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
save_fig(p2, "03_host_filter_bacterial_reads", width = 10, height = 4.5)

# ---- Figure 3: GC of bacterial contigs, kept vs removed --------------------------------------
flags <- fread(file.path(hf_dir, "contig_host_flags.tsv.gz"),
               select = c("superkingdom", "gc", "host_aligned"))[superkingdom == "Bacteria"]
flags[, status := factor(ifelse(host_aligned, "Host-aligned, removed", "Kept"),
                         levels = c("Kept", "Host-aligned, removed"))]
status_n <- table(flags$status)
status_labels <- setNames(sprintf("%s (n = %s)", names(status_n), format(as.vector(status_n), big.mark = ",")),
                          names(status_n))

p3 <- ggplot(flags, aes(gc, fill = status)) +
  geom_histogram(binwidth = 1, boundary = 0, colour = "white", size = 0.1) +
  facet_wrap(~ status, ncol = 1, scales = "free_y", labeller = as_labeller(status_labels)) +
  scale_fill_manual(values = c("Kept" = col_kept, "Host-aligned, removed" = col_removed),
                    guide = "none") +
  scale_x_continuous(limits = c(10, 90), breaks = seq(10, 90, 10)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = "GC content of Bacteria-classified contigs",
       subtitle = "Host contigs are expected near the H. vulgaris GC of ~31-33%",
       x = "Contig GC (%)", y = "Contigs") +
  theme_qc
save_fig(p3, "03_host_filter_bacterial_gc", width = 7, height = 5)

# ---- Figure 4: functions of Eukaryota-labelled ORFs on kept bacterial contigs ----------------
euk <- fread(file.path(hf_dir, "eukaryota_orfs_on_bacterial_contigs.tsv"))
if (nrow(euk)) {
  blank_na <- function(x) if (is.null(x)) NA_character_ else fifelse(x == "", NA_character_, x)
  euk[, fun := fcoalesce(blank_na(euk[["KEGGFUN"]]), blank_na(euk[["COGFUN"]]),
                         blank_na(euk[["PFAM"]]), "No functional annotation")]
  top <- euk[, .(orfs = .N, reads = sum(total_reads)), by = fun][order(-reads)]
  fwrite(top, file.path(hf_dir, "eukaryota_orfs_top_functions.tsv"), sep = "\t")

  top20 <- head(top, 20)
  top20[, fun_short := ifelse(nchar(fun) > 60, paste0(substr(fun, 1, 57), "..."), fun)]
  top20[, fun_short := factor(fun_short, levels = rev(unique(fun_short)))]

  p4 <- ggplot(top20, aes(reads / 1e3, fun_short)) +
    geom_col(width = 0.7, fill = col_kept) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(title = "Eukaryota-labelled ORFs on kept bacterial contigs",
         subtitle = sprintf("Top 20 functions by reads (all samples); %s ORFs in total, reported not removed",
                            format(nrow(euk), big.mark = ",")),
         x = "Reads (thousands)", y = NULL) +
    theme_qc + theme(panel.grid.major.y = element_blank(),
                     panel.grid.major.x = element_line(colour = "#e4e3df", size = 0.3))
  save_fig(p4, "03_host_filter_eukaryota_orf_functions", width = 8, height = 5.5)
} else {
  message("No Eukaryota-labelled ORFs on kept bacterial contigs; figure 4 skipped.")
}

# ---- Console summary -------------------------------------------------------------------------
message("\nEukaryote marker PFAM check (all contigs -> host-free bacteria):")
print(fread(file.path(hf_dir, "eukaryote_marker_pfam_check.tsv")))
message("\nPer-sample bacterial reads after filtering (lowest first):")
print(smry[order(pct_bacteria_kept),
           .(sample, day, bacteria_kept, pct_bacteria_kept = round(pct_bacteria_kept, 2),
             pct_eukaryota_orf = round(pct_eukaryota_orf, 2))])
