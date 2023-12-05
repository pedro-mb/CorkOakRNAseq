# RNA-seq analysis of cork oak Stem samples from plants after 6mo of drought (WD water deficit .vs. WW well watered)

### load raw counts table
#"RawReadCounts_stems2021.txt" available at ArrayExpress Accession: E-MTAB-13376

# Load main packages
suppressMessages(library(DESeq2))
library(ggplot2)
library(pheatmap) 
library(RColorBrewer)

countstable <- read.table("data/RawReadCounts_stems2021.txt", header=TRUE, row.names=1)

# extract gene lengths to calculate TPM counts
geneLengths <- countstable[,5]

# Remove first five columns (chr, start, end, strand, length)
countstable <- countstable[ ,6:ncol(countstable)]

# Eliminate low count genes
keep <- rowMeans(countstable) > 10

# Make a new table without the low count genes
countdata <- countstable[keep, ]

# column names from countstable:
# WW => CnC, CnP, CnX
# WD => SnC, SnP, SnX
# Phellem => CnC, SnC
# InnerBark => CnP, SnP
# Xylem => CnX, SnX
# 2, 4 and 5 are biological replicates

# Define experiment conditions (condition, tissue, etc)
coldata <- data.frame(row.names = colnames(countdata),
                      Treatment=c(rep("WW",9),rep("WD",9)),
                      Tissue=c(rep(c("Phellem", "InnerBark","Xylem"),6)))


# convert to factor
coldata$Treatment <- as.factor(coldata$Treatment)
coldata$Tissue <- as.factor(coldata$Tissue)

# order levels
coldata$Tissue <- factor(coldata$Tissue, levels = c("Phellem", "InnerBark","Xylem"))
coldata$Treatment <- relevel(coldata$Treatment, ref = "WW")

# Define colors for factor levels
ann_colors = list(Tissue = c(Phellem="#E6AB02", InnerBark="#009E73",
                             Xylem="#1F78B4"), 
                  Treatment = c(WW="#BDC3C7" ,WD="#2C3E50"))

# create deseq object
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, 
                              design= ~Tissue + Tissue:Treatment)
dds


### (OPTIONAL) Convert raw counts in Transcript per million 

# counts_to_TPM.R available here:
# https://rdrr.io/github/skimlab/CCSBUtils/src/R/counts_to_TPM.R 

source("bin/counts_to_tpm.R")
#source("counts_to_TPM.R")

# import mean read length, obtained from bam files using picardmetrics
meanReadLength_df <- read.table("data/mean_read_length_pair.txt", 
                                header = TRUE)
meanReadLength <- as.numeric(meanReadLength_df[1,])

# check if columns are in the same order
colnames(countstable)
colnames(meanReadLength_df)

GeneCounts_TPM <- counts_to_tpm(countstable, geneLengths, meanReadLength)

meanExpression_tpm <- data.frame(
  P_WW=rowMeans(GeneCounts_TPM[,c(1,4,7)]),
  P_WD=rowMeans(GeneCounts_TPM[,c(10,13,16)]),
  IB_WW=rowMeans(GeneCounts_TPM[,c(2,5,8)]),
  IB_WD=rowMeans(GeneCounts_TPM[,c(11,14,17)]),
  Xy_WW=rowMeans(GeneCounts_TPM[,c(3,6,9)]),
  Xy_WD=rowMeans(GeneCounts_TPM[,c(12,15,18)]))

# center gene expression per sample around the mean expression -> z-score
gene_means_tpm <- rowMeans(meanExpression_tpm)
gene_sd_tpm <- apply(meanExpression_tpm, 1, sd)
all_zscore_tpm <- (meanExpression_tpm - gene_means_tpm)/gene_sd_tpm


### PCAs e clustering

# Principal Component Analysis (PCA) with ggplot2

vsd <- vst(dds, blind = TRUE) 

data <- plotPCA(vsd, intgroup = c( "Treatment", "Tissue"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

pcaDATA <- plotPCA(vsd,ntop = 1000,intgroup = c("Treatment","Tissue"),
                   returnData = TRUE)

percentvar <- round(100*attr(pcaDATA,"percentVar"))

pca_plot <- ggplot(pcaDATA,aes(PC1,PC2,colour = Tissue, 
                               shape = Treatment)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ",percentvar[1],"% variance")) +
  ylab(paste0("PC2: ",percentvar[2],"% variance")) +
  scale_color_manual(values=c( "#E6AB02","#009E73", "#1F78B4")) +
  theme(axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14),
        legend.text = element_text(size = 14))

pca_plot

#ggsave(pca_plot, file = "RNAseq_Stems_PCA.png", device = "png", width = 1800, height = 1480, units = "px", dpi = 300)
#ggsave(pca_plot, file = "RNAseq_Stems_PCA.svg", width = 1800, height = 1480, units = "px", device = "svg")

sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)

rownames(sampleDistMatrix) <- paste(vsd$Treatment, vsd$Tissue, sep="-") 
colnames(sampleDistMatrix) <- colnames(vsd)

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
hdist <- pheatmap(sampleDistMatrix,
                  clustering_distance_rows=sampleDists,
                  clustering_distance_cols=sampleDists,
                  col=colors)

hdist
#ggsave(hdist, file = "RNAseq_Stems_dist.png", device = "png", width = 1800, height = 1400, units = "px", dpi = 300)
#ggsave(hdist, file = "RNAseq_Stems_dist.svg", width = 1800, height = 1400, units = "px", device = "svg")



### A: Run multigroup comparison analysis (Likelihood ratio test, ANOVA-like)

# get DEGs between tissues, i.e genes with differential expression in 
# at least one of the following pair-wise comparisons:
# Phellem vs Xylem; Phellem vs InnerBark; InnerBark vs Xylem

dds_anova <- DESeq(dds, test="LRT", reduced=~Treatment)

res <- results(dds_anova)
res <- res[order(res$padj),]

padj.cutoff <- 0.001
baseMean.cutoff <- 50
resSig <- subset(res, padj < padj.cutoff)
resSig <- subset(resSig, baseMean >= baseMean.cutoff)

# Get number of significant genes
nrow(resSig)


### Calculate mean expression of biological replicates and z-score for plotting  

# Extract Transformed Values / Regularize Logaritms 
vsd <- vst(dds,blind = FALSE)


# z-score for all genes
meanExpression <- data.frame(
  P_WW=rowMeans(assay(vsd)[,c(1,4,7)]),
  P_WD=rowMeans(assay(vsd)[,c(10,13,16)]),
  IB_WW=rowMeans(assay(vsd)[,c(2,5,8)]),
  IB_WD=rowMeans(assay(vsd)[,c(11,14,17)]),
  Xy_WW=rowMeans(assay(vsd)[,c(3,6,9)]),
  Xy_WD=rowMeans(assay(vsd)[,c(12,15,18)]))
gene_means_all <- rowMeans(meanExpression)
gene_sd_all <- apply(meanExpression, 1, sd)
all_zscore <- (meanExpression - gene_means_all)/gene_sd_all

coldata_mean <- data.frame(row.names = colnames(all_zscore),
                           Treatment=c(rep(c("WW","WD"),3)),
                           Tissue=c(rep("Phellem",2), rep("InnerBark",2),
                                    rep("Xylem", 2)))


### Clustering analysis (k-means)

source("bin/clustering_functions.R")
suppressMessages(library(tidyr))


meanExpression_cl <- data.frame(
  Phellem=rowMeans(assay(vsd)[,c(1,4,7,10,13,16)]),
  InnerBark=rowMeans(assay(vsd)[,c(2,5,8, 11,14,17)]),
  Xylem=rowMeans(assay(vsd)[,c(3,6,9, 12,15,18)]))
meanExpression_cl <- cbind(meanExpression_cl, id=rownames(assay(vsd) ))

meanExpression_DEGs <- meanExpression_cl[rownames(resSig), ] 

ncluster <- 6


meanExpression_DEGs$cluster_1 <- clusterExpr(meanExpression_DEGs[, 1:3], ncluster)$clustering


clusterDEGs_long <- meanExpression_DEGs[,c(1:5)] %>% 
  pivot_longer(c(Phellem:Xylem), names_to = "sample", values_to = "ncounts") #%>%
#mutate(ncounts = log2(ncounts) + 1)

#### block the order of x axis if necessary
levels1=c("Phellem", "InnerBark", "Xylem")

clusterDEGs_long$sample <- factor(clusterDEGs_long$sample, levels = levels1)

# Plot clusters according to expression profile
clust_1 <- ggplot(clusterDEGs_long, aes(sample, ncounts)) + 
  geom_line(aes(group = id), colour = "gray") + 
  facet_wrap(~cluster_1, scales = "free_y") +
  stat_summary(fun = "mean", geom = "line", linewidth = 1, aes(group = 1)) +
  xlab("") + ylab("Norm. read counts (vst)") +
  theme_bw() +
  theme(
    panel.grid.minor  = element_blank(), 
    panel.grid.major =  element_blank(),
    axis.text.y = element_text(size=8),
    axis.text.x = element_text(size=10),
    legend.text = element_text(size = 10),
    strip.background = element_rect(colour="black", fill=c("grey95","black"))
  )

#ggsave("DEGs_RNAseq_clust6_all.001.svg", plot = clust_1, width = 10, height = 6 )

# Plot clusters according to mean expression tendency
clust_2 <- ggplot(clusterDEGs_long, aes(sample, ncounts, group = cluster_1)) +
  facet_wrap(~cluster_1  , scales = "free_y") +
  stat_summary(fun.data = "mean_cl_boot", geom = "ribbon", fill = "cyan4", alpha = 0.5) +
  stat_summary(fun = "mean", geom = "line") + xlab("Sample") + ylab("Norm. read counts (vst)") +
  theme_bw() +
  theme(
    panel.grid.minor  = element_blank(), 
    panel.grid.major =  element_blank(),
    axis.text.y = element_text(size=10),
    axis.text.x = element_text(size=12),
    legend.text = element_text(size = 14),
    strip.background = element_rect(colour="black", fill=c("grey95","black"))
  )
#ggsave("DEGs_RNAseq_clust6_0.001.svg", plot = clust_2, width = 10, height = 6 )

#Table S3

#write.csv(merge(resSig, meanExpression_DEGs[,c(4,5)], by =0), file = "DEGs_LRT_multitissue_kmeans.csv")

### Draw heatmaps for candidate genes
mylength <- 100

cambiumMarkers <- c("LOC111988759", #WOX4 *
                    "LOC112027031", #CLE44
                    "LOC112012413", #EPFL4 *
                    "LOC112024199", #EPFL6
                    "LOC112019802", #PXY/TDR *
                    "LOC112031035", #SHR *
                    "LOC111993990", #STM.2 *
                    "LOC112032307", #STM.1 *
                    "LOC111997860", #ER  *
                    "LOC112031780", #ER?
                    "LOC111992965", #ERL1
                    "LOC112038775", #LBD4 *
                    "LOC112024503", #LBD11
                    "LOC111997151", #LBD4
                    "LOC111999314", #ARR5 *
                    "LOC112030198", #AINTEGUMENTA *
                    "LOC112015571", #KNAT/BP *
                    "LOC112038380", #CYCD3;1
                    "LOC112022222", #CYCD3;1 *
                    "LOC111985236", #AINTEGUMENTA *
                    "LOC112011953", #ARF5/MP *
                    "LOC112026395", #TMO5
                    "LOC111992348"#LHW  *
)

list <- all_zscore[rownames(all_zscore) %in% intersect(cambiumMarkers, rownames(resSig)),]

mybreaks<- c(seq(min(list), 0, length.out=ceiling(mylength/2) + 1), seq(max(list)/mylength, max(list), length.out=floor(mylength/2)))


p <- pheatmap(list, 
              cluster_rows=TRUE, 
              show_rownames=TRUE,
              cluster_cols=FALSE, 
              color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
              annotation_col = coldata_mean,
              breaks = mybreaks,
              border_color = "white",
              cellwidth = 20,        
              cellheight = 10,
              fontsize_row = 8,
              annotation_colors = ann_colors)

#ggsave("DEGsANOVA_CambiumMarkers.svg", plot = p, width = 6, height = 10.4 )


# tissue differentiation markers  
diffMarkers<- c("LOC112007487", #MYB84
                "LOC112030452", #MYB93
                "LOC111991558", #WRKY13
                "LOC112038777", #MYB83
                "LOC111996101", #KAN
                "LOC112040698", #KAN
                "LOC112010403", #ATHB15
                "LOC112031480", #PHABOLOSA Athb14
                "LOC112038274", #REVOLUTA
                "LOC112001733", #asft1
                "LOC112007140",#asft2
                "LOC112033310", #MYB93
                "LOC112028279", #GPAT5
                "LOC111990215", #APL
                "LOC112007669" #APL
)


list <- all_zscore[rownames(all_zscore) %in% intersect(diffMarkers, rownames(resSig)),]

mybreaks<- c(seq(min(list), 0, length.out=ceiling(mylength/2) + 1), seq(max(list)/mylength, max(list), length.out=floor(mylength/2)))


p<- pheatmap(list, 
             cluster_rows=TRUE, 
             show_rownames=TRUE,
             cluster_cols=FALSE, 
             color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
             annotation_col = coldata_mean,
             breaks = mybreaks,
             border_color = "white",
             cellwidth = 20,        
             cellheight = 10,
             fontsize_row = 8,
             annotation_colors = ann_colors)

#ggsave("DEGsANOVA_SpecDiffer.svg", plot = p, width = 6, height = 10.4 )



### B: Run DESeq function for pairwise DE analysis WD vs WD

dds <- DESeq(dds)

# See results names
resultsNames(dds)

res_CorkDrought <- lfcShrink(dds, coef = "TissuePhellem.TreatmentWD",type="apeglm")
res_PhloemDrought <- lfcShrink(dds, coef = "TissueInnerBark.TreatmentWD",type="apeglm")
res_XylemDrought <- lfcShrink(dds, coef = "TissueXylem.TreatmentWD",type="apeglm")

# Order by Adjusted p-value
res_CorkDrought <- res_CorkDrought[order(res_CorkDrought$padj),]
res_PhloemDrought <- res_PhloemDrought[order(res_CorkDrought$padj),]
res_XylemDrought <- res_XylemDrought[order(res_CorkDrought$padj),]

# Filter for padj < 0.05
resSig_CorkDrought <- subset(res_CorkDrought, padj < 0.05)
resSig_PhloemDrought <- subset(res_PhloemDrought, padj < 0.05)
resSig_XylemDrought <- subset(res_XylemDrought, padj < 0.05)


# Filter for |log2FC| > 1
resSig_CorkDrought <- subset(resSig_CorkDrought, abs(log2FoldChange) > 1)
resSig_PhloemDrought <-subset(resSig_PhloemDrought, abs(log2FoldChange) >1) 
resSig_XylemDrought <- subset(resSig_XylemDrought, abs(log2FoldChange) >1) 

summary(resSig_CorkDrought)
summary(resSig_PhloemDrought)
summary(resSig_XylemDrought)

# Table S4
#write.csv(resSig_CorkDrought, file="resSig_CorkDrought_padj0.05_log2FC1.csv")

# Table S5
#write.csv(resSig_PhloemDrought, file="resSig_PhloemDrought_padj0.05_log2FC1.csv")

# Table S6
#write.csv(resSig_XylemDrought, file="resSig_XylemDrought_padj0.05_log2FC1.csv")



## Get VENN diagram for DEGs between WD vs WW determined for each layer
source("bin/GOVenn_corrected.R")

resSig_cork_DvsC <- data.frame(genes=rownames(resSig_CorkDrought), logFC = resSig_CorkDrought[,2])
resSig_xylem_DvsC <- data.frame(genes=rownames(resSig_XylemDrought), logFC = resSig_XylemDrought[,2]) 
resSig_phloem_DvsC <- data.frame(genes=rownames(resSig_PhloemDrought), logFC = resSig_PhloemDrought[,2])


venn<-GOVenn(data1 = resSig_xylem_DvsC, 
             data2 = resSig_phloem_DvsC, 
             data3 = resSig_cork_DvsC , 
             plot=TRUE,
             title = "Drought DEGs", 
             label = c("Xylem",  "InnerBark","Phellem" ) , 
             circle.col = c("#66CCFF","#99FF99", "#FFCC99"))


#ggsave(plot = venn, file = "Drought_VennDEGs.svg", width = 8.2, height = 8, units = "in", device = "svg")

#ggsave("Drought_VennDEGs.svg", plot = venn, width = 6.2, height = 5.8,device = "svg")

venn_data<-GOVenn(data1 = resSig_phloem_DvsC, 
                  data2 = resSig_xylem_DvsC, 
                  data3 = resSig_cork_DvsC, 
                  plot=FALSE)

#save.image(file = "DE_workflow.RData")
#load("DE_workflow.RData")