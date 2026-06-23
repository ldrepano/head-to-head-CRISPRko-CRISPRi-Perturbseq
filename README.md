
This repository presents the analysis of analogous CRISPRko and CRISPRi Perturb-seq screens conducted by the Genetic Perturbation Platform and Genomics Platform at the Broad Institute of MIT and Harvard. 

To set up this analysis, use the following commands in terminal

```
python3 -m venv ~/head-to-head-CRISPRko-CRISPRi-perturbseq
source ~/head-to-head-CRISPRko-CRISPRi-perturbseq/bin/activate
pip install -r requirements.txt
```

The analyze-singlecelldata directory includes all R and WDL scripts employed to analyze the data at single-cell resolution (in contrast to the SCEPTRE results, which summarize perturbation effects per guide). This data is not included in this repository due to size constraints. Raw single-cell RNA-seq data as fastq files have been deposited at [GEO](link pending). Aligned, quality-filtered Perturb-seq datasets with guide assignments are available as RDS files to [Zenodo](https://doi.org/10.5281/zenodo.20722064). 