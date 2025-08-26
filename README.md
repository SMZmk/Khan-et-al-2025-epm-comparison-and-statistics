# arabis-drought-tolerance-epm-analyses

This repository contains the R and Python scripts used to analyze gene expression and regulatory networks related to drought tolerance in the plant species **Arabis nemorensis** and **Arabis sagittata**. The analyses leverage deep learning models to identify **Expression Predictive Motifs (EPMs)** and investigate their distribution and enrichment within differentially expressed genes.

---

## Overview

The code in this repository corresponds to the analyses described in the paper titled "Polygenic basis of molecular mechanisms underlying drought tolerance differences in ecologically different grassland species."

The core purpose of these scripts is to:

1.  **Identify and Cluster EPMs**: Group similar EPMs into clusters to represent functional regulatory units.
2.  **Enrichment Analysis**: Determine which EPM clusters are significantly enriched in genes that are up- or down-regulated in response to drought stress.
3.  **Functional Annotation**: Use Gene Ontology (GO) enrichment to annotate the biological functions associated with the enriched EPM clusters.
4.  **Interspecies Comparison**: Compare the EPMs and their associated regulatory networks between the two species to explain their different drought tolerance strategies.

---

## Repository Contents

* **`EPM_DGEenrichment_analyses.v3.2.R`**: This script performs the primary enrichment analysis. It calculates various statistical metrics (Z-scores, Poisson tests) to determine if **Expression Predictive Motifs (EPMs)** are enriched in differentially expressed genes. The output includes tables and plots used for the main figures and supplementary materials.
* **`EPM_JASPAR_venn-analyses.v1.R`**: This script identifies transcription factor binding sites (**TFBS**) by comparing EPMs against the JASPAR database. It then uses Venn diagrams (Euler plots) to visualize the overlap and uniqueness of these TFBS across different species and stress conditions.
* **`EPM_Quadrant_enrichment_2025F.v2.1.R`**: A script for comparative enrichment analysis. It tests whether EPM clusters are enriched in specific gene expression quadrants (e.g., genes that are highly expressed in one species but not the other) to identify species-specific regulatory patterns.
* **`GO-EPMs_species_comp.R`**: This script focuses on associating GO terms with EPM clusters to provide functional context. It aggregates and compares GO enrichment results between the two species to highlight differences in the biological processes activated during drought.
* **`epm_to_reference_alingment.0.py`**: A Python script designed to create multi-sequence FASTA alignments. It takes a reference alignment, a set of pattern sequences, and a mapping file to generate a combined, padded alignment file. This is useful for visualizing sequence motifs in a conserved context.
