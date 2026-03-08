# DEG Analysis of GSE66597

## Overview
This repository contains a differential gene expression (DEG) analysis of the **GSE66597** dataset from the Gene Expression Omnibus (GEO) The analysis compares gene expression profiles between **control U251 astrocyte cells** and **cells infected with Influenza A H5N1 virus**.

The aim of this project is to identify genes genes that show significant expression changes due to H5N1 infection and to understand the cellular response at the transcriptomic level.

---

## Dataset Information
- **Dataset ID:** GSE66597  
- **Organism:** *Homo sapiens*  
- **Cell line:** U251 astrocyte cell line  
- **Platform:** GPL570 (Affymetrix Human Genome U133 Plus 2.0 Array)  
- **Total Samples:** 18

### Experimental Design

**Control group**
- con6h (3 replicates)
- con12h (3 replicates)
- con24h (3 replicates)

**H5N1 infection group**
- hm6h (3 replicates)
- hm12h (3 replicates)
- hm24h (3 replicates)

Each group represents gene expression measurements at **6, 12, and 24 hours post infection (hpi)**.

---

## Analysis Method

Differential gene expression analysis was performed using **limma** through the GEO2R platform.

Analysis parameters:

- Statistical method: **limma**
- Multiple testing correction: **Benjamini–Hochberg (FDR)**
- Significance criteria:
  - Adjusted p-value (adj.P.Val) < 0.05
  - |log2 Fold Change| > 1

---

## Results

The analysis identified **4545 differentially expressed genes (DEGs)** between control and H5N1-infected astrocyte cells.

Visualization results include:
- Volcano Plot
- Mean Difference Plot
- DEG Summary Diagram

These findings suggest that **H5N1 infection induces extensive transcriptomic changes in human astrocytes**, which may be associated with inflammatory responses, cellular stress, and neurological complications.

---

## Repository Structure
- Dataset Gse66597_raw.txt
- README.md
- Report.pdf
- Result Table GSE66597.txt
- Script GSE66597.R

** File Description**
- `Dataset GSE66597_RAW.txt` → Raw dataset used for analysis  
- `Result Table GSE66597.txt` → Differential gene expression results  
- `Script GSE66597.R` → R script used for the analysis  
- `Report.pdf` → Final analysis report  
- `README.md` → Project documentation  

---

## Tools and Software

- R Programming Language
- GEOquery
- limma
- GEO2R (NCBI GEO)

---

## References

NCBI Gene Expression Omnibus (GEO)  
https://www.ncbi.nlm.nih.gov/geo/

Dataset GSE66597  
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66597
