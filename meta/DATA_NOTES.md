# DATA_NOTES — SqueezeMeta output, noninduced HVBMT samples

Survey date: 2026-09-29. Source: `/mnt/hdd3/sqm_noninduced/`. **This directory is read-only**; it was inspected only with `du`, `ls`, `head` and header comparisons.

## TL;DR
- All 30 sample IDs in the SqueezeMeta tables match `meta/HVBMT_metadata_noninduced.csv` exactly, **in the same order**. No renaming is needed.
- **The data is mostly host.** In every sample, *Hydra vulgaris* is the most abundant taxon at every rank. Only **1.4–9.4 % of reads are Bacteria**.
- The KO / COG / PFAM tables in `results/tables/` **include host ORFs**, so bacterial functional analysis must first subset to Bacteria (see [Host dominance](#host-dominance--the-main-caveat)).
- For analysis, use the files in `results/tables/` (≤ 9 MB each, apart from the per-ORF and per-contig taxonomy files). Do not load `data/`, `intermediate/` or `temp/`.

## Run info
- SqueezeMeta **v1.8.0**, project name `sqm_noninduced`, coassembly mode. The run started 2026-09-22 and finished 2026-09-24. All steps 1–21 and `sqm2tables.py` completed (`progress` ends with `END`).
- Assembly: MEGAHIT, contigs < 200 bp removed. 1,145,849 contigs, 825 Mb in total, N50 = 824 bp, longest contig 269,752 bp.
- ORFs: Prodigal. RNAs: Barrnap and Aragorn. Homology search: Diamond against GenBank nr, eggNOG (COG) and KEGG; HMMER3 against Pfam.
- Taxonomy: LCA. 29.5 % of contigs are classified at superkingdom level and 14.9 % at genus level.
- Mapping: Bowtie2. Binning: MetaBAT2, then DAS Tool, with CheckM2 for bin quality. Pathways: MinPath (KEGG, MetaCyc).
- The full methods text is in `methods.txt`; parameters are in `parameters.pl` and `SqueezeMeta_conf.pl`.

## Directory map (`du -h --max-depth`)
```
653G  /mnt/hdd3/sqm_noninduced/
638G  ├── data/
551G  │   ├── raw_fastq/          rRNA-filtered reads (*.nonrrna.fastq)
 65G  │   ├── bam/                one BAM + .bai per sample (sqm_noninduced.HVBMT##.bam)
 23G  │   ├── megahit/            assembly
      │   └── 00.sqm_noninduced.samples   sample → fastq map (60 lines = 30 samples × 2 pairs)
6.6G  ├── intermediate/         binners 969M, checkm 69M, checkm2 33M, minpath
4.4G  ├── results/
851M  │   ├── tables/             ← SUMMARY TABLES FOR ANALYSIS (48 files)
121M  │   ├── bins/
      │   └── 01–21.* files       per-step outputs (see below)
2.8G  ├── temp/
3.4M  ├── ext_tables/           12.*.cog.stamp, 12.*.kegg.stamp (STAMP-format)
738M  ├── sqm_noninduced.zip
      └── creator.txt, methods.txt, parameters.pl, SqueezeMeta_conf.pl, progress, syslog
```

### `results/` top level: large files that should not be loaded casually
| File | Size | Notes |
|---|---|---|
| `01.sqm_noninduced.fasta` | 865 MB | contigs |
| `13.sqm_noninduced.orftable` | 738 MB | per-ORF annotation + per-sample counts. `SQMtools::loadSQM()` reads this |
| `19.sqm_noninduced.contigtable` | 703 MB | per-contig table |
| `03.sqm_noninduced.{gff,fna,faa}` | 504 / 416 / 228 MB | ORF predictions |
| `06.*.fun3.tax.wranks` | 61–70 MB | ORF taxonomy |
| `12.*.{cog,kegg}.funcover` | 36–43 MB | |

Small `results/` files that are useful: `10.sqm_noninduced.mappingstat` (1.5 KB), `21.sqm_noninduced.stats` (7 KB), `18.sqm_noninduced.bintable` (5 KB), `20.sqm_noninduced.{kegg,metacyc}.pathways`, `02.sqm_noninduced.16S.txt`.

## Summary tables (`results/tables/`, prefix `sqm_noninduced.`)
Each table has one row per feature and one column per sample (30 samples). The first column is the feature ID and has no header.

**Taxonomy abundance.** Files are named `{superkingdom,phylum,class,order,family,genus,species}.{allfilter,nofilter,prokfilter}.abund.tsv`.
- Values are raw read counts. Row names are full lineage strings, e.g. `k_Bacteria;p_Pseudomonadota;...`.
- Row counts: phylum.prokfilter has 144 taxa and genus.prokfilter has 1,930.
- `superkingdom.*` also includes `k_No CDS`, `k_Unclassified` and `k_Unmapped` rows. Their totals therefore equal all reads.

**Function.** Files are named `{KO,COG,PFAM}.{abund,tpm,cov,bases,copyNumber}.tsv`.
- Feature counts: KO 11,762, COG 21,336, PFAM 9,199.
- `KO.names.tsv` and `COG.names.tsv` hold the columns `Name` and `Path` (the pathway hierarchy, `|`-separated). PFAM row IDs already include the name, e.g. `PF00001 [7 transmembrane receptor ...]`.

**Per-ORF and per-contig taxonomy (large).** `orf.tax.{all,no,prok}filter.tsv` (~128–130 MB) and `contig.tax.*.tsv` (~134–136 MB).

**Other.** `bin.tax.tsv`, `orf.16S.tsv`, `RecA.tsv` and `orf.marker.genes.tsv` (30 MB).

### First lines of the key tables (columns truncated)
`phylum.prokfilter.abund.tsv`
```
                                        HVBMT01 HVBMT21 HVBMT41 HVBMT04 HVBMT24 HVBMT44 ...
k_Archaea;p_Candidatus Thermoplasmatota     210     244     344     383     267     318 ...
k_Archaea;p_Candidatus Thorarchaeota          0       0       0       2       0       0 ...
k_Archaea;p_Candidatus Woesearchaeota         0       0       0       0       0       0 ...
```
`genus.prokfilter.abund.tsv`
```
k_Archaea;p_Candidatus Thermoplasmatota;c_Candidatus Poseidoniia;...;g_Unclassified Candidatus Poseidoniia  204 244 344 371 ...
k_Archaea;p_Candidatus Thermoplasmatota;c_Unclassified ...;g_Unclassified Candidatus Thermoplasmatota       6   0   0  12 ...
```
`KO.abund.tsv` (raw counts) and `KO.tpm.tsv`
```
        HVBMT01 HVBMT21 HVBMT41 HVBMT04 ...          HVBMT01  HVBMT21  HVBMT41 ...
K00001      358     505     877     269 ...   K00001  3.659    5.587    8.807  ...
K00002     1492    1220     907    1517 ...   K00002 27.225   26.328   19.188  ...
K00003      112     107     172      72 ...   K00003  1.315    1.568    2.028  ...
```
`KO.names.tsv`
```
        Name                                        Path
K00001  alcohol dehydrogenase [EC:1.1.1.1]          Metabolism; Carbohydrate metabolism; Glycolysis / Gluconeogenesis | ...
K00002  alcohol dehydrogenase (NADP+) [EC:1.1.1.2]  Metabolism; Carbohydrate metabolism; Glycolysis / Gluconeogenesis | ...
```
`COG.abund.tsv`
```
COG0001   441   515   758   300 ...
COG0002   152   205   334   102 ...
```
`PFAM.abund.tsv`: note that the top rows are eukaryotic GPCR domains, i.e. host.
```
PF00001 [7 transmembrane receptor (rhodopsin family)]    8221   6948   6443 ...
PF00002 [7 transmembrane receptor (Secretin family)]    12934  10922   9860 ...
```

## Sample-name check against metadata: PASS
- The metadata file `meta/HVBMT_metadata_noninduced.csv` has 30 rows plus a header. Columns are `SampleID, Replicate, Temp, Time, Treatment, GroupCombined, sample`. There are no duplicate IDs and no CRLF line endings.
- The design matches CLAUDE.md: Day 0 at 18 °C × 3 replicates, plus 3 temperatures (8, 18, 28 °C) × Days 1/7/28 × 3 replicates, for 30 samples.
- **All 36 per-sample tables** have exactly those 30 IDs as columns, **in the same order as the metadata rows**. These are all 24 taxonomy `abund` tables, plus KO/COG/PFAM × `abund`, `tpm`, `cov`, `bases` and `copyNumber`. Every taxonomy table has an identical header.
- The same 30 IDs appear in `data/00.sqm_noninduced.samples`, `data/bam/` and `10.mappingstat`.
- ⚠️ `21.stats` lists samples in sorted ID order, not metadata order. Always join on sample name, never on column position.
- Replicates are coded by the tens digit: replicate 1 = 0x/1x, replicate 2 = 2x/3x, replicate 3 = 4x/5x. For example, 1_28 is HVBMT02/22/42.

## Sequencing depth and mapping (`10.mappingstat`)
- Reads per sample range from 37.4 M (HVBMT56) to 81.7 M (HVBMT12); the total is 1.63 G reads.
- Mapping rate to the coassembly is 96.3–99.5 % for all samples except **HVBMT32 (8 °C, D7, replicate 2): 91.4 %**. It also has the most unmapped reads (8.6 %).
- Depth varies about 2-fold between samples, so counts must be normalized before any comparison.

## Host dominance — the main caveat
Per-sample read fractions, from `superkingdom.nofilter.abund.tsv`:

| Sample | Temp/Day | Bacteria % | Eukaryota % | Unclassified % | Unmapped % |
|---|---|---|---|---|---|
| HVBMT01 | 18/D0 | 3.6 | 51.9 | 42.3 | 1.7 |
| HVBMT21 | 18/D0 | 4.4 | 51.9 | 41.4 | 1.7 |
| HVBMT41 | 18/D0 | 4.3 | 50.6 | 40.8 | 3.7 |
| HVBMT04 | 18/D1 | 2.0 | 46.5 | 49.3 | 1.9 |
| HVBMT24 | 18/D1 | 3.4 | 45.4 | 50.1 | 0.9 |
| HVBMT44 | 18/D1 | 3.8 | 44.9 | 50.1 | 0.9 |
| HVBMT02 | 28/D1 | 2.8 | 41.7 | 54.1 | 1.2 |
| HVBMT22 | 28/D1 | 8.0 | 51.2 | 39.7 | 0.8 |
| HVBMT42 | 28/D1 | 4.0 | 46.8 | 47.7 | 1.2 |
| HVBMT06 | 8/D1 | 2.1 | 45.2 | 50.8 | 1.7 |
| HVBMT26 | 8/D1 | 3.8 | 49.1 | 45.8 | 1.0 |
| HVBMT46 | 8/D1 | 4.6 | 52.2 | 41.7 | 1.1 |
| HVBMT10 | 18/D7 | 5.7 | 52.0 | 41.4 | 0.5 |
| HVBMT30 | 18/D7 | 7.3 | 50.9 | 40.5 | 0.7 |
| HVBMT50 | 18/D7 | 4.1 | 48.8 | 44.8 | 1.0 |
| HVBMT08 | 28/D7 | 2.7 | 53.0 | 43.4 | 0.7 |
| HVBMT28 | 28/D7 | 3.4 | 51.0 | 43.9 | 1.1 |
| HVBMT48 | 28/D7 | 3.2 | 48.3 | 47.0 | 1.1 |
| HVBMT12 | 8/D7 | 3.0 | 42.9 | 50.5 | 2.5 |
| HVBMT32 | 8/D7 | 6.0 | 40.5 | 43.6 | 8.6 |
| HVBMT52 | 8/D7 | 9.4 | 47.5 | 40.1 | 1.3 |
| HVBMT16 | 18/D28 | 2.7 | 50.2 | 45.5 | 1.3 |
| HVBMT36 | 18/D28 | 2.3 | 47.2 | 49.3 | 0.9 |
| HVBMT56 | 18/D28 | 2.6 | 39.4 | 56.8 | 0.8 |
| HVBMT14 | 28/D28 | 2.7 | 48.9 | 46.9 | 1.1 |
| HVBMT34 | 28/D28 | 3.3 | 52.8 | 42.5 | 0.9 |
| HVBMT54 | 28/D28 | 4.0 | 49.2 | 45.1 | 1.3 |
| HVBMT18 | 8/D28 | 1.5 | 51.7 | 45.0 | 0.8 |
| HVBMT38 | 8/D28 | 1.4 | 50.1 | 46.7 | 0.9 |
| HVBMT58 | 8/D28 | 2.0 | 51.6 | 43.6 | 1.2 |

Archaea and Viruses each make up less than 0.2 %. The "No CDS" category is omitted from the table.

**What this means for analysis**
1. According to `21.stats`, the most abundant taxon in every sample is *Hydra vulgaris* (Cnidaria) at every rank. The Unclassified fraction (40–57 %) is probably mostly host as well, e.g. poorly annotated Hydra transcripts.
2. The **function tables (`KO/COG/PFAM.*.tsv`) are community-wide and include host ORFs.** This is evident from the GPCR domains at the top of the PFAM table. They **cannot be used as-is** for any of the bacterial questions (objectives 1–5).
3. Recommended route: in R, `SQMtools::loadSQM("/mnt/hdd3/sqm_noninduced")`, then `subsetTax(SQM, "superkingdom", "Bacteria")`, and use the functional tables of that bacterial subset. `loadSQM()` reads the 738 MB orftable and 703 MB contigtable, so check memory first. Loading once and saving the bacterial subset to an `.rds` is worthwhile. `loadSQMlite()` on `results/tables/` **cannot** do this subsetting, because it has no per-ORF taxonomy.
4. The bacterial fraction varies about 7-fold, from 1.4 % in HVBMT38 to 9.4 % in HVBMT52. Bacterial functional profiles must therefore be normalized within the bacterial subset, e.g. TPM or relative abundance recomputed after subsetting, or per-taxon normalization. They should not be normalized to total reads.
5. The lowest bacterial fractions are in **8 °C Day 28 (1.4–2.0 %)**. This could itself be a biological signal (objective 1C), but it also means those samples have the fewest bacterial reads and therefore the noisiest bacterial profiles.
6. The `prokfilter` taxonomy tables are the relevant ones for community composition. Even so, the first rows of the phylum and genus tables are Archaea, so filter explicitly to `k_Bacteria` if the focus is bacteria only.

**Dominant bacterial phyla** (reads summed across all samples, nofilter): Spirochaetota 24.9 M > Pseudomonadota 21.5 M >> unclassified Bacteria 9.0 M > Bacteroidota 2.0 M > Actinomycetota 1.7 M > Cyanobacteriota 0.65 M > Myxococcota 0.31 M.

## Bins (`18.bintable`, DAS Tool, 7 bins; 3 are ≥ 50 % complete)
| Bin | Taxonomy | Size | GC % |
|---|---|---|---|
| metabat2.1 | Bacteria; Pseudomonadota; Gammaproteobacteria; Pseudomonadaceae | 4.5 Mb | 64.3 |
| metabat2.6 | Bacteria; Spirochaetota; *Leptospira* | 3.8 Mb | 40.7 |
| metabat2.5 | Bacteria; Alphaproteobacteria; *Phyllobacterium* | 2.8 Mb | 59.7 |
| metabat2.4 | *Hydra vulgaris* | 5.0 Mb | 31.4 |
| metabat2.2 | Eukaryota (unclassified; 16S says Bacteria) | 104.7 Mb | 28.0 |
| metabat2.3 | Eukaryota | 0.4 Mb | 41.5 |
| metabat2.7 | Unclassified | 0.4 Mb | 42.0 |

The three bacterial bins (Pseudomonadaceae, *Leptospira* and *Phyllobacterium*) match the two dominant phyla. They are candidates for genome-resolved expression analysis. Completeness and contamination values are in `18.bintable`.

## Project-side state
- `meta/scripts/` and `meta/result/figures/` do not exist yet, and `meta/result/` is empty. The first script (`01_load.R`) should create them, following the CLAUDE.md conventions.
- The other metadata files in `meta/` (`HVBMT_metadata_induced.csv`, and `HVBMT_metadata_combined_updated.csv`, which covers both datasets) are not used for this SqueezeMeta run.

## Reproduce (read-only commands used)
```bash
D=/mnt/hdd3/sqm_noninduced; T=$D/results/tables/sqm_noninduced
du -h --max-depth=1 $D; du -h --max-depth=1 $D/data $D/intermediate $D/results
ls -la $D/results $D/results/tables
head -4 $T.KO.abund.tsv | cut -c1-400
cat $D/results/10.sqm_noninduced.mappingstat
# sample-name check: header of every per-sample table vs metadata SampleID column
tail -n +2 meta/HVBMT_metadata_noninduced.csv | cut -d, -f1 > meta_ids.txt
for f in $D/results/tables/*.{abund,tpm,cov,bases,copyNumber}.tsv; do
  head -1 $f | tr '\t' '\n' | grep '^HVBMT' | diff -q - meta_ids.txt >/dev/null && echo "$f OK" || echo "$f DIFF"
done
# per-sample superkingdom %
awk -F'\t' 'NR==1{for(i=2;i<=NF;i++)s[i]=$i;n=NF;next}{for(i=2;i<=n;i++){v[$1,i]=$i;t[i]+=$i}}
  END{for(i=2;i<=n;i++)printf "%s\t%.1f\t%.1f\n",s[i],100*v["k_Bacteria",i]/t[i],100*v["k_Eukaryota",i]/t[i]}' \
  $T.superkingdom.nofilter.abund.tsv
```
