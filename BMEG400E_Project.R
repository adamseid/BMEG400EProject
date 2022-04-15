### INSTALL PACKAGES ###
BiocManager::install("affyio")
BiocManager::install("affxparser")
BiocManager::install("GEOquery")
install.packages("GEOquery")
BiocManager::install("affy")
BiocManager::install("affxparser")
BiocManager::install("DESeq2")
install.packages("DESeq2")

### LOAD LIBRARIES ###
library(affyio)
library(affxparser)
library(GEOquery)
library(affy)
library(stats)
library(data.table)
library(ggbiplot)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library( "DESeq2" )

### Unzip g-files ###
gunzip("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434456_0IU_CD4_2.rma-gene-default.chp.gz", remove=FALSE)
gunzip("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434457_0IU_CD4_3.rma-gene-default.chp.gz", remove=FALSE)
gunzip("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434458_0IU_CD4_4.rma-gene-default.chp.gz", remove=FALSE)
gunzip("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434459_0IU_CD4_5.rma-gene-default.chp.gz", remove=FALSE)
gunzip("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434460_0IU_CD4_6.rma-gene-default.chp.gz", remove=FALSE)

gunzip("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434461_10IU_CD4_7.rma-gene-default.chp.gz", remove=FALSE)
gunzip("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434462_10IU_CD4_8.rma-gene-default.chp.gz", remove=FALSE)
gunzip("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434463_10IU_CD4_9.rma-gene-default.chp.gz", remove=FALSE)
gunzip("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434464_10IU_CD4_10.rma-gene-default.chp.gz", remove=FALSE)
gunzip("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434465_10IU_CD4_11.rma-gene-default.chp.gz", remove=FALSE)
gunzip("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434466_10IU_CD4_12.rma-gene-default.chp.gz", remove=FALSE)

gunzip("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434467_NF_CD4_13.rma-gene-default.chp.gz", remove=FALSE)
gunzip("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434468_NF_CD4_14.rma-gene-default.chp.gz", remove=FALSE)
gunzip("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434469_NF_CD4_15.rma-gene-default.chp.gz", remove=FALSE)
gunzip("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434470_NF_CD4_16.rma-gene-default.chp.gz", remove=FALSE)
gunzip("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434471_NF_CD4_17.rma-gene-default.chp.gz", remove=FALSE)
gunzip("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434472_NF_CD4_18.rma-gene-default.chp.gz", remove=FALSE)

### Load Chp files###
OIUReplicate2 = readChp("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434456_0IU_CD4_2.rma-gene-default.chp",withQuant = TRUE)
OIUReplicate3 = readChp("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434457_0IU_CD4_3.rma-gene-default.chp",withQuant = TRUE)
OIUReplicate4 = readChp("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434458_0IU_CD4_4.rma-gene-default.chp",withQuant = TRUE)
OIUReplicate5 = readChp("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434459_0IU_CD4_5.rma-gene-default.chp",withQuant = TRUE)
OIUReplicate6 = readChp("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434460_0IU_CD4_6.rma-gene-default.chp",withQuant = TRUE)

IOIUReplicate1 = readChp("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434461_10IU_CD4_7.rma-gene-default.chp",withQuant = TRUE)
IOIUReplicate2 = readChp("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434462_10IU_CD4_8.rma-gene-default.chp",withQuant = TRUE)
IOIUReplicate3 = readChp("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434463_10IU_CD4_9.rma-gene-default.chp",withQuant = TRUE)
IOIUReplicate4 = readChp("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434464_10IU_CD4_10.rma-gene-default.chp",withQuant = TRUE)
IOIUReplicate5 = readChp("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434465_10IU_CD4_11.rma-gene-default.chp",withQuant = TRUE)
IOIUReplicate6 = readChp("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434466_10IU_CD4_12.rma-gene-default.chp",withQuant = TRUE)

NFReplicate1 = readChp("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434467_NF_CD4_13.rma-gene-default.chp",withQuant = TRUE)
NFReplicate2 = readChp("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434468_NF_CD4_14.rma-gene-default.chp",withQuant = TRUE)
NFReplicate3 = readChp("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434469_NF_CD4_15.rma-gene-default.chp",withQuant = TRUE)
NFReplicate4 = readChp("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434470_NF_CD4_16.rma-gene-default.chp",withQuant = TRUE)
NFReplicate5 = readChp("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434471_NF_CD4_17.rma-gene-default.chp",withQuant = TRUE)
NFReplicate6 = readChp("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434472_NF_CD4_18.rma-gene-default.chp",withQuant = TRUE)

### Get all OIU IDs ###
OIUProb = c(as.numeric(OIUReplicate2$QuantificationEntries$ProbeSetName),
            as.numeric(OIUReplicate3$QuantificationEntries$ProbeSetName),
            as.numeric(OIUReplicate4$QuantificationEntries$ProbeSetName),
            as.numeric(OIUReplicate5$QuantificationEntries$ProbeSetName),
            as.numeric(OIUReplicate6$QuantificationEntries$ProbeSetName))

### Get all OIU Values ###
OIUVal = c(OIUReplicate2$QuantificationEntries$QuantificationValue,
           OIUReplicate3$QuantificationEntries$QuantificationValue,
           OIUReplicate4$QuantificationEntries$QuantificationValue,
           OIUReplicate5$QuantificationEntries$QuantificationValue,
           OIUReplicate6$QuantificationEntries$QuantificationValue)

### Get all 1OIU Values ###
IOIUVal = c(IOIUReplicate1$QuantificationEntries$QuantificationValue,
            IOIUReplicate2$QuantificationEntries$QuantificationValue,
            IOIUReplicate3$QuantificationEntries$QuantificationValue,
            IOIUReplicate4$QuantificationEntries$QuantificationValue,
            IOIUReplicate5$QuantificationEntries$QuantificationValue,
            IOIUReplicate6$QuantificationEntries$QuantificationValue)

### Get all 1OIU Values ###
NFVal = c(NFReplicate1$QuantificationEntries$QuantificationValue,
          NFReplicate2$QuantificationEntries$QuantificationValue,
          NFReplicate3$QuantificationEntries$QuantificationValue,
          NFReplicate4$QuantificationEntries$QuantificationValue,
          NFReplicate5$QuantificationEntries$QuantificationValue,
          NFReplicate6$QuantificationEntries$QuantificationValue)

### Get list of all data ###
fullData <- data.frame("0IU2" = OIUReplicate2$QuantificationEntries$QuantificationValue, "0IU3" = OIUReplicate3$QuantificationEntries$QuantificationValue, 
                       "0IU4" = OIUReplicate4$QuantificationEntries$QuantificationValue, "0IU5" = OIUReplicate5$QuantificationEntries$QuantificationValue, 
                       "0IU6" = OIUReplicate6$QuantificationEntries$QuantificationValue, "10IU1" = IOIUReplicate1$QuantificationEntries$QuantificationValue, 
                       "10IU2" = IOIUReplicate2$QuantificationEntries$QuantificationValue, "10IU3" = IOIUReplicate3$QuantificationEntries$QuantificationValue, 
                       "10IU4" = IOIUReplicate4$QuantificationEntries$QuantificationValue, "10IU5" = IOIUReplicate5$QuantificationEntries$QuantificationValue, 
                       "10IU6" = IOIUReplicate6$QuantificationEntries$QuantificationValue, "NF1" = NFReplicate1$QuantificationEntries$QuantificationValue,
                       "NF2" = NFReplicate2$QuantificationEntries$QuantificationValue, "NF3" = NFReplicate3$QuantificationEntries$QuantificationValue,
                       "NF4" = NFReplicate4$QuantificationEntries$QuantificationValue, "NF5" = NFReplicate5$QuantificationEntries$QuantificationValue,
                       "NF6" = NFReplicate6$QuantificationEntries$QuantificationValue)

### Get Names of each Study ###
names <- c("0IU2", "0IU3", "0IU4", "0IU5", "0IU6", "10IU1","10IU2", "10IU3", "10IU4", "10IU5", 
           "10IU6", "NF1","NF2", "NF3","NF4", "NF5","NF6")
fullDataTwo.studies <- c(rep("0IU", 5), rep("10IU",6), rep("2IU", 6))

### Transpose Matrix ###
fullDataTwo <- transpose(fullData)
row.names(fullDataTwo) <- names
colnames(fullDataTwo) <- IOIUReplicate1$QuantificationEntries$ProbeSetName

### Remove un-needed variables ###
rm(OIUReplicate2,OIUReplicate3,OIUReplicate4,OIUReplicate5,OIUReplicate6,IOIUReplicate1,
   IOIUReplicate2,IOIUReplicate3,IOIUReplicate4, IOIUReplicate5, IOIUReplicate6,
   NFReplicate1,NFReplicate2,NFReplicate3,NFReplicate4,NFReplicate5,NFReplicate6, OIUProb,
   fullData.pca, res.pca, fullDataTwo.studies, IOIUVal, NFVal, OIUVal)

### Perform PCA Analysis ###
# center is a flag which ensures the variables within "fullDataTwo" is shifted to 0
# Scale flag ensures the variables within "fullDataTwo" is scaled to contain variance.
fullDataTwo.pca <- prcomp(fullDataTwo, center = TRUE,scale = TRUE)
# ellipse flag draws a normal data ellipse for the 0IU, 10IU and 2IU groups
# var.axes was set to false to remove arraws within our PCA plot. This is to clean up the graph
# as the arrows were so numerous it covered the entire plot (including the pca points)
# labels were set to be each individual study for 0IU, 10IU, and 2IU. There are 17 studies in total
# 5 done on 0IU, 6 done on 10IU and 6 done on 2IU. Thus we will have 17 labels in total
# groups was used to group the 17 studies into 3 overall groups (0IU, 10IU, 2IU). This was done
# to better represent the spread of the PCA plot and how seperated each cluster was
ggbiplot(fullDataTwo.pca,ellipse=TRUE,var.axes=FALSE,labels=rownames(fullDataTwo), groups=fullDataTwo.studies)

### Create HeatMap ###
# This function transforms our data into a matrix. The reason this was done was because our 
# heatmap accepts only matrix
fullData.matrix <- data.matrix(fullData)
row.names(fullData.matrix) <- row.names(fullData)
colnames(fullData.matrix) <- colnames(fullData)
# This function center's our data within fullData.matrix. This is important as it 
# improves the accuracy of our heatmap since it normalizes our data
fullData.matrix.scaled = t(apply(fullData.matrix, 1, scale))
# This function creates our heat map. Unlike the normal heat map function offered by r, 
# this function calculates a complex heat map. This is important as it allows us to create
# a heat map based on the expression of our data (this is what the author of the paper did).The
# col tag determines the colors that will be used on our heat map. We chose a gradient from red
# to blue to mimic the literature heat map. We made tags show_column_names and show_row_dend false.
# These tags are responsible of showing column names and row dendrogram. We made it false because
# they required a lot of computational power to produce and did not affect the overall graph (they
# were not needed).
Heatmap(fullData.matrix.scaled, name = "expression", row_km = 5,
        col = colorRamp2(c(-2, 2), c("red", "blue")),
        show_column_names = FALSE, row_title = NULL, show_row_dend = FALSE)

### Create Differenetial Expression ###
ratData <- read.delim("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/RaGeneVsRat230_BestMatch.txt")
head(ratData)
fullDataThree <- data.matrix(fullDataTwo)
names <- colnames(fullDataTwo)
for (x in 1:nrow(fullDataThree)) {
  for(y in 1:ncol(fullDataThree)){
    fullDataThree[x,y] = as.integer(fullDataThree[x,y])
  }
}
dds <- DESeqDataSetFromMatrix(fullDataThree, DataFrame(names), ~ 1)
dds <- DESeq(dds)
res <- results(dds, contrast=c("0IU", "10IU", "2IU"))[order(res$pvalue),]
resSig <- subset(res, padj < 0.06)
# Find 5 downregulated data
head(resSig[ order( resSig$log2FoldChange ), ])
# Find 5 upregulated data
head(resSig[ order( resSig$log2FoldChange, decreasing=TRUE), ])

### BASH SCRIPT ###
### Install Library ###
conda install -c bioconda/label/cf201901 apt-probeset-summarize
### Create chp files from CEL files ###
# Perform a rma analysis because this is what the author used. Each rma analysis has a 
# log2 transformation perfomed on it. The output-dir specifies the cel we will be using. The
# --cc-chp-output creates the chp file
apt-probeset-summarize -a rma -d chip.cdf -o --cc-chp-output output-dir GSM2434456_0IU_CD4_2.CEL
apt-probeset-summarize -a rma -d chip.cdf -o --cc-chp-output output-dir GSM2434457_0IU_CD4_3.CEL
apt-probeset-summarize -a rma -d chip.cdf -o --cc-chp-output output-dir GSM2434458_0IU_CD4_4.CEL
apt-probeset-summarize -a rma -d chip.cdf -o --cc-chp-output output-dir GSM2434459_0IU_CD4_5.CEL
apt-probeset-summarize -a rma -d chip.cdf -o --cc-chp-output output-dir GSM2434460_0IU_CD4_6.CEL
apt-probeset-summarize -a rma -d chip.cdf -o --cc-chp-output output-dir GSM2434461_10IU_CD4_7.CEL
apt-probeset-summarize -a rma -d chip.cdf -o --cc-chp-output output-dir GSM2434462_10IU_CD4_8.CEL
apt-probeset-summarize -a rma -d chip.cdf -o --cc-chp-output output-dir GSM2434463_10IU_CD4_9.CEL
apt-probeset-summarize -a rma -d chip.cdf -o --cc-chp-output output-dir GSM2434464_10IU_CD4_10.CEL
apt-probeset-summarize -a rma -d chip.cdf -o --cc-chp-output output-dir GSM2434465_10IU_CD4_11.CEL
apt-probeset-summarize -a rma -d chip.cdf -o --cc-chp-output output-dir GSM2434466_10IU_CD4_12.CEL
apt-probeset-summarize -a rma -d chip.cdf -o --cc-chp-output output-dir GSM2434467_NF_CD4_13.CEL
apt-probeset-summarize -a rma -d chip.cdf -o --cc-chp-output output-dir GSM2434468_NF_CD4_14.CEL
apt-probeset-summarize -a rma -d chip.cdf -o --cc-chp-output output-dir GSM2434469_NF_CD4_15.CEL
apt-probeset-summarize -a rma -d chip.cdf -o --cc-chp-output output-dir GSM2434470_NF_CD4_16.CEL
apt-probeset-summarize -a rma -d chip.cdf -o --cc-chp-output output-dir GSM2434471_NF_CD4_17.CEL
apt-probeset-summarize -a rma -d chip.cdf -o --cc-chp-output output-dir GSM2434472_NF_CD4_18.CEL
