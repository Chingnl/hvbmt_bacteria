# HVBMT bacterial metatranscriptome project: workflow

Update the "Current state" section and the log at the end of every session.

## Current state
- **Last updated:** 2026-09-30
- **Stage:** Host-filter scripts 01–03 written and syntax-checked; not yet run.
- **Next steps:** host-filter plan (details below)
  - [ ] Run `bash meta/scripts/01_contigs_vs_hydra.sh` (minimap2, about 32 threads)
  - [ ] Run `Rscript meta/scripts/02_host_filter.R` (about 5 GB RAM)
  - [ ] Run `Rscript meta/scripts/03_host_filter_qc.R` and review the figures and `host_filter_summary.txt`
  - [ ] Record the results in DATA_NOTES.md and in this file
- **Open questions / blockers:**
  - The host-alignment thresholds (≥ 50 % coverage, ≥ 90 % identity) may need tuning after the QC.
  - Once the Eukaryota-ORF report exists, decide whether to drop those ORFs.

### Host-filter plan (2026-09-30)
**Why:** KneadData had already removed host reads with bowtie2 (strain 105 index, `~/db/kneaddata/hvul105-bt2-db/`). Even so, Eukaryota still make up about 40–53 % of reads and Bacteria only 1.4–9.4 %. bowtie2 is not splice-aware, so re-mapping the reads to the same index would mostly repeat that work. The filter therefore acts on the SqueezeMeta contigs instead. Nothing is written under `/mnt/hdd3/`.

**Step 1: `01_contigs_vs_hydra.sh`**
- Input: `/mnt/hdd3/sqm_noninduced/results/01.sqm_noninduced.fasta` (all 1.15 M contigs) and `~/db/kneaddata/hvul105-bt2-db/HydraT2T_105_genomic.fna`.
- Run `minimap2 -x splice:hq` (spliced preset for accurate sequences) (conda env `minimap2`, about 32 threads).
- Output: `meta/result/host_filter/contigs_vs_hydra105.paf`.
- Aligning all contigs, not just the Bacteria ones, also shows how much of the 550 Mb of unclassified contigs is host.

**Step 2: `02_host_filter.R`**
1. `SQM <- readRDS("/mnt/hdd3/sqm_data_cache.rds")` (the full `loadSQM()` object, about 4.8 GB of RAM).
2. Flag a contig as host if its alignments to Hydra cover ≥ 50 % of its length at ≥ 90 % identity. 90 % is enough because the cultures are strain 105, the same strain as the reference.
3. `subsetTax(SQM, "superkingdom", "Bacteria")`, then `subsetContigs()` to remove the host-flagged contigs. Rescale TPM within the final subset.
4. **Report, do not drop:** ORFs whose own ORF-level taxonomy is `k_Eukaryota` but that sit on the remaining bacterial contigs. Write a table with the ORF ID, contig ID, ORF and contig taxonomy, KEGG/COG/PFAM annotation, GC and per-sample read counts.
5. Also report, without dropping, any bacterial contigs that sit in the host bins metabat2.4 (*Hydra*) and metabat2.2 (Eukaryota, 104.7 Mb).
6. Outputs in `meta/result/host_filter/`:
   - `sqm_bacteria_hostfree.rds`;
   - a per-contig flag table;
   - the Eukaryota-ORF report;
   - per-sample bacterial reads before and after removing host contigs, plus the reads held in Eukaryota-labelled ORFs.

**Step 3: `03_host_filter_qc.R`** (figures to `meta/result/figures/` as PDF, PNG and SVG)
- Per-sample stacked bars: Bacteria kept, Bacteria removed as host, Eukaryota, Unclassified.
- GC histogram of bacterial contigs kept vs removed. Removed contigs are expected to sit near the host GC of about 31–33 %.
- Sanity check: eukaryotic PFAM domains (for example the PF00001/PF00002 GPCRs) should be about zero after filtering. Any left over are tracked in the Eukaryota-ORF report.
- Summary of the Eukaryota-ORF report: number of ORFs, share of bacterial reads, and top functions.
- Bacterial read depth per sample, especially 8 °C Day 28 (1.4–2.0 % of reads are Bacteria).

## Key decisions
| Date | Decision | Rationale |
|---|---|---|
| 2026-09-29 | Subset to Bacteria before any functional analysis | Host dominates reads |
| 2026-09-30 | Host-filter the contigs with splice-aware minimap2, not the reads | Reads were already bowtie2-filtered against Hydra 105 by KneadData; bowtie2 misses spliced host transcripts |
| 2026-09-30 | Host-alignment identity threshold ≥ 90 % | Cultures are *H. vulgaris* strain 105, the same strain as the reference genome (not AEP) |
| 2026-09-30 | Report, do not drop, Eukaryota-labelled ORFs on Bacteria-labelled contigs | Inspect them first; they may be host contamination or horizontally transferred genes |
| 2026-09-30 | Load `/mnt/hdd3/sqm_data_cache.rds` instead of running `loadSQM()` | The cache is the full SQM object (about 4.8 GB of RAM) |

## Task log (newest first)

### 2026-09-30: Host-filter scripts written
- **Output:**
  - `meta/scripts/01_contigs_vs_hydra.sh`: minimap2 `-x splice:hq --secondary=no` → `meta/result/host_filter/contigs_vs_hydra105.paf`.
  - `meta/scripts/02_host_filter.R`: `sqm_bacteria_hostfree.rds`, per-contig flags, read accounting, Eukaryota-ORF and host-bin reports, PFAM marker check.
  - `meta/scripts/03_host_filter_qc.R`: four figures (`03_host_filter_*`) plus a top-functions table.
- **Status:** syntax-checked only; not run.
- **Notes:**
  - `svglite` is not installed, so the SVGs are written with base `grDevices::svg`.
  - Eukaryote marker PFAMs checked: PF00001/2/3 (GPCRs), PF00125 (histone), PF00022 (actin).

### 2026-09-30: Host-filter planning and notes update
- **Output:** host-filter plan (above); `CLAUDE.md` and `meta/DATA_NOTES.md` updated. DATA_NOTES.md was trimmed to 183 lines by summarizing the example table rows.
- **Findings:**
  - KneadData had already filtered reads against the strain 105 bowtie2 index, yet host still dominates.
  - `sqm_data_cache.rds` is the full SQM object.
  - Contig taxonomy: Bacteria 91.7 k contigs / 68 Mb / mean GC 53 %; Eukaryota 245 k / 206 Mb / 33 %; Unclassified 808 k / 550 Mb / 31 %.
  - About 17 k Bacteria contigs have 20–39 % GC, which overlaps the host range.

### 2026-09-29: Data survey
- **Output:** `meta/DATA_NOTES.md`
- **Findings:** host fraction ~40–53%, bacteria ~1.4–9.4% (source of rules 1–2)
