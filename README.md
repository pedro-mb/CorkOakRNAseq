# CorkOakRNAseq

This repository contains all the detailed commands used for differential gene expression analysis, based on RNA-seq from cork oak stem tissues grown under water deficit (WD) and well watered (WW) conditions. 

**Publication title:** Drought impact on phellem development: identification of novel gene regulators and evidence of photosynthetic activity
**DOI:** https://doi.org/10.1101/2023.12.26.573371

Raw RNA-seq data and read count matrix are available at ArrayExpress (Accession: E-MTAB-13376)
Methods:
RNA samples for each tissue (Phellem, InnerBark and Xylem), collected from stems grown under WW and WD were submitted to high-throughput sequencing (100 bp, paired-end reads) using DNBSEQTM at BGI Tech. (Hong Kong). Read quality was evaluated using FastQC (v0.11.9). Reads were pre-processed to remove adapters and low-quality bases using Trimmomatic v0.39 (Bolger et al., 2014) with paired-end mode and additional parameters: LEADING:3, TRAILING:3, SLIDINGWINDOW:4:15, MINLEN:36. High-quality reads were mapped against Q. suber reference genome GCF_002906115.1 (CorkOak1.0) (Ramos et al., 2018), using STAR v2.7.7a (Dobin et al., 2013) with modified parameters (twopassMode = Basic, alignIntronMax = 100 000). Read counts per gene were obtained with featureCounts (Liao et al., 2014), using exon as feature. 

Differential expression analysis described here was conducted using the DESeq2 (Love et al., 2014) in R v.4.2.1. Genes with mean read count < 10 for all samples were removed. Principal component analysis (PCA) was conducted using normalized counts of the 1000 genes showing the highest variance. Differential expression analysis between WD and WW conditions and the three tissues was performed using a multifactorial design (~Tissue + Tissue: Treatment) and the ‘apeglm’ log fold-change shrinkage estimator (Zhu et al., 2019). Differentially expressed genes (DEGs) found for WD vs. WW were filtered for adjusted P-value (Padj) ≤ 0.05 and fold-change (|log2FC| ≥ 1). DE analysis between the three tissues was conducted using the Likelihood Ratio Test and DEGs were determined by adjusted P-value (Padj) ≤ 0.001 and baseMean > 50. DEGs between tissues were grouped using k-means clustering. Candidate A. thaliana orthologues were identified using Blastp analysis of cork oak amino acid sequences on TAIR10 (Araport11) predicted proteins.

