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

### Unzip g-files ###
gunzip("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434456_0IU_CD4_2.rma-gene-default.chp.gz", remove=FALSE)
gunzip("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434457_0IU_CD4_3.rma-gene-default.chp.gz", remove=FALSE)
gunzip("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434458_0IU_CD4_4.rma-gene-default.chp.gz", remove=FALSE)
gunzip("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434459_0IU_CD4_5.rma-gene-default.chp.gz", remove=FALSE)
gunzip("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434460_0IU_CD4_6.rma-gene-default.chp.gz", remove=FALSE)

gunzip("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434456_0IU_CD4_2.CEL.gz", remove=FALSE)
gunzip("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434457_0IU_CD4_3.CEL.gz", remove=FALSE)
gunzip("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434458_0IU_CD4_4.CEL.gz", remove=FALSE)
gunzip("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434459_0IU_CD4_5.CEL.gz", remove=FALSE)
gunzip("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434460_0IU_CD4_6.CEL.gz", remove=FALSE)

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

### Load Chp files ###
OIUReplicate2 = readChp("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434456_0IU_CD4_2.rma-gene-default.chp",withQuant = TRUE)
OIUReplicate3 = readChp("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434457_0IU_CD4_3.rma-gene-default.chp",withQuant = TRUE)
OIUReplicate4 = readChp("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434458_0IU_CD4_4.rma-gene-default.chp",withQuant = TRUE)
OIUReplicate5 = readChp("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434459_0IU_CD4_5.rma-gene-default.chp",withQuant = TRUE)
OIUReplicate6 = readChp("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434460_0IU_CD4_6.rma-gene-default.chp",withQuant = TRUE)

OIUReplicate2Data = read.celfile("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434456_0IU_CD4_2.CEL",intensity.means.only=FALSE)
OIUReplicate3Data = read.celfile("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434457_0IU_CD4_3.CEL",intensity.means.only=FALSE)
OIUReplicate4Data = read.celfile("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434458_0IU_CD4_4.CEL",intensity.means.only=FALSE)
OIUReplicate5Data = read.celfile("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434459_0IU_CD4_5.CEL",intensity.means.only=FALSE)
OIUReplicate6Data = read.celfile("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434460_0IU_CD4_6.CEL",intensity.means.only=FALSE)

length <- length(OIUReplicate2Data$OUTLIERS) + length(OIUReplicate3Data$OUTLIERS) + 
  length(OIUReplicate4Data$OUTLIERS) + length(OIUReplicate5Data$OUTLIERS) + length(OIUReplicate6Data$OUTLIERS)
OIUReplicate2Data$MASKS
a <- read.probematrix("/Users/ethio/OneDrive/Desktop/UBC Spring 2022/BMEG 400E/Assignments/project/DataFiles/GSM2434456_0IU_CD4_2.CEL")
OIUReplicate3Data[[1]]
OIUReplicate3Data

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
length(OIUProb)
### Get all OIU Values ###
OIUVal = c(OIUReplicate2$QuantificationEntries$QuantificationValue,
           OIUReplicate3$QuantificationEntries$QuantificationValue,
           OIUReplicate4$QuantificationEntries$QuantificationValue,
           OIUReplicate5$QuantificationEntries$QuantificationValue,
           OIUReplicate6$QuantificationEntries$QuantificationValue, l)
l <- rep(0, 29214)
InternationalUnit <-rep(c(0),each=146070)
OIUdf <- data.frame(ID = OIUProb, Value = OIUVal, Type = InternationalUnit)

### Get all 1OIU IDs ###
IOIUProb = c(as.numeric(IOIUReplicate1$QuantificationEntries$ProbeSetName),
             as.numeric(IOIUReplicate2$QuantificationEntries$ProbeSetName),
             as.numeric(IOIUReplicate3$QuantificationEntries$ProbeSetName),
             as.numeric(IOIUReplicate4$QuantificationEntries$ProbeSetName),
             as.numeric(IOIUReplicate5$QuantificationEntries$ProbeSetName),
             as.numeric(IOIUReplicate6$QuantificationEntries$ProbeSetName))

### Get all 1OIU Values ###
IOIUVal = c(IOIUReplicate1$QuantificationEntries$QuantificationValue,
            IOIUReplicate2$QuantificationEntries$QuantificationValue,
            IOIUReplicate3$QuantificationEntries$QuantificationValue,
            IOIUReplicate4$QuantificationEntries$QuantificationValue,
            IOIUReplicate5$QuantificationEntries$QuantificationValue,
            IOIUReplicate6$QuantificationEntries$QuantificationValue)

InternationalUnit <-rep(c(1),each=175284)
IOIUdf <- data.frame(ID = IOIUProb, Value = IOIUVal, Type = InternationalUnit)

### Get all NF IDs ###
NFProb = c(as.numeric(NFReplicate1$QuantificationEntries$ProbeSetName),
           as.numeric(NFReplicate2$QuantificationEntries$ProbeSetName),
           as.numeric(NFReplicate3$QuantificationEntries$ProbeSetName),
           as.numeric(NFReplicate4$QuantificationEntries$ProbeSetName),
           as.numeric(NFReplicate5$QuantificationEntries$ProbeSetName),
           as.numeric(NFReplicate6$QuantificationEntries$ProbeSetName))

### Get all 1OIU Values ###
NFVal = c(NFReplicate1$QuantificationEntries$QuantificationValue,
          NFReplicate2$QuantificationEntries$QuantificationValue,
          NFReplicate3$QuantificationEntries$QuantificationValue,
          NFReplicate4$QuantificationEntries$QuantificationValue,
          NFReplicate5$QuantificationEntries$QuantificationValue,
          NFReplicate6$QuantificationEntries$QuantificationValue)

InternationalUnit <-rep(c(2),each=175284)
NFdf <- data.frame(ID = NFProb, Value = NFVal, Type = InternationalUnit)

### Get list of all data ###
fullData <- rbind(OIUdf,IOIUdf)
fullData <- rbind(fullData,NFdf)
fullData <- data.frame("0IU2" = OIUReplicate2$QuantificationEntries$QuantificationValue, "0IU3" = OIUReplicate3$QuantificationEntries$QuantificationValue, 
                       "0IU4" = OIUReplicate4$QuantificationEntries$QuantificationValue, "0IU5" = OIUReplicate5$QuantificationEntries$QuantificationValue, 
                       "0IU6" = OIUReplicate6$QuantificationEntries$QuantificationValue, "10IU1" = IOIUReplicate1$QuantificationEntries$QuantificationValue, 
                       "10IU2" = IOIUReplicate2$QuantificationEntries$QuantificationValue, "10IU3" = IOIUReplicate3$QuantificationEntries$QuantificationValue, 
                       "10IU4" = IOIUReplicate4$QuantificationEntries$QuantificationValue, "10IU5" = IOIUReplicate5$QuantificationEntries$QuantificationValue, 
                       "10IU6" = IOIUReplicate6$QuantificationEntries$QuantificationValue, "NF1" = NFReplicate1$QuantificationEntries$QuantificationValue,
                       "NF2" = NFReplicate2$QuantificationEntries$QuantificationValue, "NF3" = NFReplicate3$QuantificationEntries$QuantificationValue,
                       "NF4" = NFReplicate4$QuantificationEntries$QuantificationValue, "NF5" = NFReplicate5$QuantificationEntries$QuantificationValue,
                       "NF6" = NFReplicate6$QuantificationEntries$QuantificationValue)
fullData <- data.frame("0IU" = OIUVal, "10IU" = IOIUVal, 
                       "NF" = NFVal)

fullData.pca <- prcomp(fullData, center = TRUE,scale. = TRUE)
summary(fullData.pca)
biplot(prcomp(fullData, scale = TRUE))
plot(fullData.pca)
