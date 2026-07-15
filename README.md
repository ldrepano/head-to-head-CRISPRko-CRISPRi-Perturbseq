
This repository hosts all code required for the analysis presented in "Direct comparison of CRISPR knockout and interference with Perturb-seq
" [(pre-print)](https://www.biorxiv.org/content/10.64898/2026.07.04.736492v1) (Broad Institute of MIT and Harvard). 

To set up this analysis, use the following commands in terminal

```
python3 -m venv ~/head-to-head-CRISPRko-CRISPRi-perturbseq
source ~/head-to-head-CRISPRko-CRISPRi-perturbseq/bin/activate
pip install -r requirements.txt
```

The analyze-singlecelldata directory includes all R and WDL scripts employed to analyze the data at single-cell resolution. The analyze-SCEPTRE-results directory (Python) uses data that represents perturbation effects (relative to cells with negative control guides) aggregated across all cells carrying the same guide per sample. 

Aligned, quality-filtered Perturb-seq datasets with guide assignments have been deposited as RDS files to Zenodo at doi.org/10.5281/zenodo.2072206446. The results of differential expression analysis (SCEPTRE) per guide for all samples have been deposited to Zenodo at doi.org/10.5281/zenodo.21136747
