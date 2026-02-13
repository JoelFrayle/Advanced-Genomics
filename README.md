# Advanced Genomics — Project 5: Replication-driven copy number variation detectable in RNA-seq

**Authors:** Joel Frayle Moreno, Francesco Mazza  
**Date:** 13-02-2026  

## Interactive HTML report (GitHub Pages)
GitHub’s file viewer cannot render large HTML files, so the report is published with GitHub Pages:

➡️ **Live report:** https://joelfrayle.github.io/Advanced-Genomics/

## Repository contents
- `index.html` — rendered interactive report (published on GitHub Pages).
- `Report.rmd` — source R Markdown used to generate the report.

## Project overview
This project tests the hypothesis that **DNA replication induces copy-number variation along the bacterial chromosome**, and that this gradient can be detected in **RNA-seq gene expression**.

### Main pipeline
1. Load a gene-by-sample expression matrix in **logTPM** (`log_tpm.csv`) and sample metadata (`sample_table.csv`).
2. Download *E. coli* genome annotation (**GFF**) from NCBI and extract gene coordinates.
3. Map expression rows to genomic positions via gene identifiers.
4. Bin genes into fixed **100 kb windows** and map bin midpoints onto circular coordinates \(\theta \in [0,2\pi)\).
5. Compute mean logTPM per bin per sample.
6. Fit a per-sample circular linear model:
   \[
   y(\theta) = \alpha + \beta \cos(\theta) + \gamma \sin(\theta)
   \]
   optionally weighting by the number of genes per bin.
7. Derive effect size and direction:
   - Amplitude \(A=\sqrt{\beta^2+\gamma^2}\)
   - Phase \(\phi=\mathrm{atan2}(\gamma,\beta)\)
   - Naive contrast \(\Delta = 2A\)
8. Identify **oriC** using **dnaA** as a proxy and compute \(\theta_{ori}\).
9. Test peak directionality (Rayleigh + V-tests) and define an ori-anchored contrast:
   \[
   \Delta_{ori} = y(\theta_{ori}) - y(\theta_{ori}+\pi)
   = 2(\beta\cos\theta_{ori} + \gamma\sin\theta_{ori})
   \]
10. Integrate measured growth rates and evaluate association + prediction (correlations + repeated 5-fold cross-validation).

## Notes on reproducibility
To reproduce the report locally, the dataset files used by the analysis must be present alongside `Report.rmd`, including (at minimum):
- `log_tpm.csv`
- `sample_table.csv`

Main R packages used include: `dplyr`, `data.table`, `stringr`, `scales`, `tibble`, plus `circular`, `CircStats`, `ggplot2`, and `pROC` in later sections.
