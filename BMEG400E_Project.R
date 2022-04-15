### INSTALL PACKAGES ###
BiocManager::install("affyio")
BiocManager::install("affxparser")
BiocManager::install("GEOquery")
install.packages("GEOquery")
BiocManager::install("affy")

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
   fullData.pca, res.pca, names, fullDataTwo.studies, IOIUVal, NFVal, OIUVal)

### Perform PCA Analysis ###
fullDataTwo.pca <- prcomp(fullDataTwo, center = TRUE,scale. = TRUE)
ggbiplot(fullDataTwo.pca,ellipse=TRUE,var.axes=FALSE,labels=rownames(fullDataTwo), groups=fullDataTwo.studies)

### Create HeatMap ###
fullData.matrix <- data.matrix(fullData)
row.names(fullData.matrix) <- row.names(fullData)
colnames(fullData.matrix) <- colnames(fullData)
base_mean = rowMeans(fullData.matrix)
fullData.matrix.scaled = t(apply(fullData.matrix, 1, scale))
type = gsub("s\\d+_", "", colnames(fullData.matrix))
ha = HeatmapAnnotation(type = type, annotation_name_side = "left")
Heatmap(fullData.matrix.scaled, name = "expression", row_km = 5,
        col = colorRamp2(c(-2, 2), c("red", "blue")),
        show_column_names = FALSE, row_title = NULL, show_row_dend = FALSE)
