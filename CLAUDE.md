# HVBMT bacterial metatranscriptome project

## Data (READ-ONLY — never modify, move, or delete anything here)
- SqueezeMeta output: /mnt/hdd3/sqm_noninduced/  (653 GB; mostly intermediate files)
- Use summary tables in results/ and results/tables/ for analysis
- Metadata: ~/HVBMT_project/meta/HVBMT_metadata_noninduced.csv

## Rules for working with large files
- Never cat/print whole data files. Use head, wc -l, du -sh, or read the first rows in R/Python.
- Check file sizes before loading anything.

## Conventions
- Language: R with SQMtools + ggplot2
- Environment: conda env base
- Scripts go in ~/HVBMT_project/meta/scripts/, numbered (01_load.R, 02_taxa_barplot.R, ...)
- Figures go in ~/HVBMT_project/meta/result/figures/ as both PDF, PNG, and SVG
- Every figure comes from a script. No one-off console plotting.

## Biology context
- Samples: labeled as HVBMT## (Hydra Vulgaris Bacterial MetaTranscriptome),
- Sample description: bacterial community cDNA purified from H. vulgaris 
- Sample conditions: 18C (control), 28C (heat stress), 8C (cold stress), Noninduced (not Mitomycin C-treated)
- Sample timepoints: collected repeatedly from the same culture at Day 1, 7 and 28 with 3 biological replicates for each temperature; Day 0 only has 18C (3 replicates)
- Sample total: 30 samples
- Key objectives: 
    1) Acute vs. Chronic Thermal Adaptation Dynamics
        1A) What is the trajectory of the initial shock response (Day 1) versus long-term acclimation (Day 28)?
        1B) At 28°C and 8°C, do transcriptional stress responses spike acutely at Day 1 and resolve by Day 28, or does chronic exposure sustain expression of alternative metabolic/survival pathways?
        1C) Is there community-level stabilization or functional collapse over 28 days?
        1D) Do the transcriptomes converge toward an alternate stable state by Day 28, or does long-term thermal divergence persist?
    2) Biofilm, Motility, and Host–Microbe Attachment Mechanics
        2A) Is bacterial sessility favored at control conditions and disrupted under thermal stress?
        2B) How do environmental sensing systems coordinate the response?
    3) Oxidative, DNA, and Protein Stress Response Mechanisms
        3A) How do the detoxification, DNA damage/repair, and chaperone networks partition between heat and cold?( follow up question: which pathways are transcribed under 28°C versus 8°C?)
        3B) Are translation and RNA modification actively prioritized under cold stress? (follow up question: How do aminoacyl-tRNA synthetases and ribonuclease enzymes respond transcriptionally at 8°C to maintain protein synthesis efficiency?)
    4) Metabolic Shift in Nutrition Sources
        4A) Are general carbon, nitrogen, sulfur, and lipid metabolic pathways stably expressed across all temperature and timepoints?
        4B) Do metatrascriptomes of heat or cold stressed communities prioritize one metabolic pathway over another?
    5) Virulence and Host Immune Evasion
        5A) Does temperature affect expression of virulence genes or other genes that aid in evading host immunity?    