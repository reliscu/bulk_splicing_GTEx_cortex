library(edgeR)
library(data.table)

source("/mnt/lareaulab/reliscu/code/SampleNetwork/SampleNetwork_1.08.r")

setwd("/mnt/lareaulab/reliscu/projects/NSF_GRFP/analyses/bulk/GTEx/cortex")

datExprT <- fread("GTEx_cortex_counts_All_501_outliers_removed_TMM_ComBat_SMGEBTCH_corrected.csv", data.table=FALSE)

sampleinfo1 <- fread("/mnt/lareaulab/reliscu/projects/NSF_GRFP/data/bulk/GTEx/cortex/GTEx_cortex_sampleinfo.csv", data.table=FALSE)
sampleinfo1$Mean_age <- sapply(strsplit(sampleinfo1$AGE, "-"), function(x) mean(as.numeric(x)))
sampleinfo1$grouplabels1 <- "All"

# Order columns by variable they will be colored by:

covar <- "SMTSD" # Tissue type
sampleinfo1 <- sampleinfo1[order(sampleinfo1[,covar]),]
datExprT <- datExprT[, c(1, match(sampleinfo1[,1], colnames(datExprT)))]

# Normalize counts

tmmf <- calcNormFactors(datExprT[,-1], method="TMM")
datExprT[,-1] <- sweep(datExprT[,-1], MARGIN=2, FUN="/", STATS=tmmf)

skip1 <- 1 # An integer describing the number of feature information columns
indices1 <- list(seq(2, ncol(datExprT)))
samplelabels1 <- 1 # An integer that points to the column number in sampleinfo1 containing the sample labels that will appear in plots.  Note: these 					sample labels must be identical to the sample column headers in datExprT or an error will be triggered.
grouplabels1 <- grep("grouplabels1", colnames(sampleinfo1))
subgroup1 <- grep(covar, colnames(sampleinfo1))

projectname1 <- "GTEx_cortex_counts_TMMF"

btrait1 <- c(
  grep("SMCENTER", colnames(sampleinfo1)), # BSS collection site
  grep("SMNABTCH", colnames(sampleinfo1)), # Batch when DNA/RNA was isolated and extracted from a sample
  grep("SMGEBTCH", colnames(sampleinfo1)), # Batch when DNA/RNA from a sample was analyzed
  grep("SEX", colnames(sampleinfo1)),
  grep("SMTSD", colnames(sampleinfo1)),
  grep("DTHHRDY", colnames(sampleinfo1)),
  grep("SMRIN", colnames(sampleinfo1)),
  grep("SMMAPRT", colnames(sampleinfo1)), # Mapping Rate: Proportion of Mapped Reads to Total Reads
  grep("SMEXNCRT", colnames(sampleinfo1)), # High Quality Exonic Rate: Proportion of Exonic Reads among High Quality Reads
  grep("SMRRNART", colnames(sampleinfo1)), # rRNA Rate: Proportion of rRNA Reads among Mapped Reads
  grep("SMDPMPRT", colnames(sampleinfo1)), # Duplicate Rate of Mapped: Proportion of Mapped Duplicate Reads to Mapped Reads
  grep("Mean_age", colnames(sampleinfo1))
)

# SMNABTCHT: The process by which DNA/RNA was isolated
# SMGEBTCHT: The data on which DNA/RNA was analyzed

asfactors1 <- c(
  grep("SMCENTER", colnames(sampleinfo1)),
  grep("SMNABTCH", colnames(sampleinfo1)),
  grep("SMGEBTCH", colnames(sampleinfo1)),
  grep("SEX", colnames(sampleinfo1)),
  grep("SMTSD", colnames(sampleinfo1)),
  grep("DTHHRDY", colnames(sampleinfo1))
)

SampleNetwork(
  datExprT=datExprT,
  method1="correlation",
  impute1=FALSE,
  subset1=NULL,
  skip1=skip1,
  indices1=indices1,
  sampleinfo1=sampleinfo1,
  subgroup1=subgroup1,
  samplelabels1=samplelabels1,
  grouplabels1=grouplabels1,
  fitmodels1=TRUE,
  whichmodel1="univariate",
  whichfit1="pc1",
  btrait1=btrait1,
  trait1=NULL,
  asfactors1=NULL,
  projectname1=projectname1,
  cexlabels1=0.7,
  normalize1=FALSE,
  replacenegs1=FALSE,
  exportfigures1=TRUE,
  verbose=TRUE
)
