#!/usr/bin/env bash
# 01_contigs_vs_hydra.sh
# Align every SqueezeMeta contig to the Hydra vulgaris strain 105 T2T genome with a
# splice-aware aligner, so host transcripts that slipped past KneadData's (unspliced)
# bowtie2 filter can be flagged in 02_host_filter.R.
#
# Input (read-only):
#   /mnt/hdd3/sqm_noninduced/results/01.sqm_noninduced.fasta   (1.15 M contigs, 865 MB)
#   ~/db/kneaddata/hvul105-bt2-db/HydraT2T_105_genomic.fna     (same genome KneadData used)
# Output:
#   meta/result/host_filter/contigs_vs_hydra105.paf
#   meta/result/host_filter/minimap2.log
#
# Usage: bash meta/scripts/01_contigs_vs_hydra.sh      (THREADS=16 to override; FORCE=1 to overwrite)

set -euo pipefail

CONTIGS=/mnt/hdd3/sqm_noninduced/results/01.sqm_noninduced.fasta
HYDRA=$HOME/db/kneaddata/hvul105-bt2-db/HydraT2T_105_genomic.fna
OUT_DIR=$HOME/HVBMT_project/meta/result/host_filter
PAF=$OUT_DIR/contigs_vs_hydra105.paf
LOG=$OUT_DIR/minimap2.log
THREADS=${THREADS:-32}

for f in "$CONTIGS" "$HYDRA"; do
  [[ -s $f ]] || { echo "Missing input: $f" >&2; exit 1; }
done
if [[ -s $PAF && ${FORCE:-0} != 1 ]]; then
  echo "$PAF already exists; set FORCE=1 to overwrite." >&2
  exit 1
fi
mkdir -p "$OUT_DIR"

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate minimap2

{
  echo "date:     $(date -Is)"
  echo "minimap2: $(minimap2 --version)"
  echo "contigs:  $CONTIGS"
  echo "genome:   $HYDRA"
  echo "threads:  $THREADS"
} > "$LOG"

# -x splice:hq    spliced alignment preset for accurate sequences (assembled transcripts)
# --secondary=no  one placement per contig region is enough to flag host
# Write to a temp file first so an interrupted run never leaves a partial PAF behind.
minimap2 -x splice:hq --secondary=no -t "$THREADS" "$HYDRA" "$CONTIGS" \
  > "$PAF.tmp" 2>> "$LOG"
mv "$PAF.tmp" "$PAF"

echo "alignments:      $(wc -l < "$PAF")" >> "$LOG"
echo "contigs aligned: $(cut -f1 "$PAF" | sort -u | wc -l)" >> "$LOG"
echo "done:            $(date -Is)" >> "$LOG"
tail -3 "$LOG"
