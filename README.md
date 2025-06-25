# Circadian-rhythm-GWAS-by-Circ-deviation-score
This repository contains scripts and workflows of Dissect the genetic architecture of circadian rhythms using a quantitative circadian deviation score

---

## Table of Contents

- [Overview](#overview)  
- [Methods](#methods)  
  - [Circadian Gene Assembly](#circadian-gene-assembly)  
  - [GTEx Genotype Preprocessing and Deviation Score Calculation](#gtex-genotype-preprocessing-and-deviation-score-calculation)  
  - [SNP Mapping and Association Testing](#snp-mapping-and-association-testing)  
  - [Confounding Factor Control](#confounding-factor-control)  
- [Results](#results)  
  - [Circadian Variation in Gene Expression](#circadian-variation-in-gene-expression)  
  - [Genome-wide Mapping of Circ-SNPs](#genome-wide-mapping-of-circ-snps)  
  - [Functional Characterization of Circ-SNPs](#functional-characterization-of-circ-snps)  
  - [Circadian Trait Associations](#circadian-trait-associations)  
  - [Drug Repurposing](#drug-repurposing)  
- [License](#license)  

---

## Overview

This project identifies SNPs associated with circadian gene expression deviations using weighted scores derived from baboon circadian gene data and human GTEx datasets. Downstream analyses explore their functional roles, genomic distribution, connection with known circadian traits, and potential therapeutic implications.

---

## Methods

### Circadian Gene Assembly

- **Script**: `explore_baboon.R`  
  Extracts circadian genes from olive baboon data and assigns p-value-based weights to compute deviation scores.

### GTEx Genotype Preprocessing and Deviation Score Calculation

- **Scripts**:  
  - `remove_head.sh`, `partition_gt.sh`, `cat.sh`: Divide GTEx genotype data into manageable chunks.  
  - `GTEx_extract_GT_auto.R`, `GTEx_extract_GT_sex.R`: Impute missing genotypes on autosomes and sex chromosomes (run in parallel).  
  - `cleanup_gt.R`, `cleanup_gt_sex.R`: Clean imputed genotype data.  
  - `granulate_gt.R`: Further split processed genotype files to reduce runtime.  
  - `GTEx_deviation_by_tissue_weighted_complete.R`: Normalize gene expression, compute Z-scores, and calculate circadian deviation scores.

### SNP Mapping and Association Testing

- **Scripts**:  
  - `regression.R`: Run parallelized linear regressions for association testing.  
  - `combine_regression_result_weighted_complete.R`: Collect significant SNPs.  
  - `compile_sumstats.R`: Compile summary statistics per tissue.

### Confounding Factor Control

- **Script**: `genotype_corr.R`  
  Validates that tissue sampling time does not confound association signals.

---

## Results

### Circadian Variation in Gene Expression

- **Script**: `dev_score.R`  
  Visualizes circadian deviation scores across all tissues.

### Genome-wide Mapping of Circ-SNPs

- **Scripts**:  
  - `GTEx_significant_SNPs_info.R`: Annotates significant SNPs, checks for overlap with circadian gene sets and GTEx QTLs.  
  - `sumstats_plot.R`: Generates Q-Q and Manhattan plots.  
  - `genomic_inflation_factor.R`: Calculates genomic inflation factors.  
  - `circ_snp_stat.R`: Summarizes the number of associations, sample sizes, and circadian gene counts per tissue.

### Functional Characterization of Circ-SNPs

- **Scripts**:  
  - `minor_allele_dev_score.R`: Analyzes deviation scores for minor alleles vs. controls.  
  - `circ_snp_ancestral_allele.R`: Retrieves ancestral alleles and shows that major Circ-SNP alleles tend to be ancestral.

- **Distribution on Genome**:  
  - `manhattan.R`: Visualizes genome-wide distribution including chrX.  
  - `tissue_overlap_upset.R`: Plots overlap of significant SNPs across tissues.  
  - `snp_location.R`, `snp_location2.R`: Show genomic location and 3' transcript enrichment of Circ-SNPs.

### Circadian Trait Associations

- **Scripts**:  
  - `GWAS_comparison_snp.R`: Cross-references Circ-SNPs with SNPs linked to circadian traits.  
  - `GWAS_comparison_gene.R`: Compares Circ-SNP-harboring genes with circadian-trait genes and highlights Circ-regulators.

### Drug Repurposing

- **Script**: `drug_repurposing.R`  
  Analyzes therapeutic areas targeted by drugs interacting with Circ-regulators and compares them with all DrugBank entries. Enrichment analysis highlights drug categories related to circadian-disrupted conditions.

---

## License

This project is for academic and research use. Please cite the corresponding publication when using the code.  
License terms will be added upon publication.

