library(ChIPQC)

setwd("/mnt/iacchus/joe/raw_data/Thymus/Fezf2_ChIP_seq/")

samples <- read.csv("sample_sheet.csv")
chip_obj <- ChIPQC(samples, annotation = "mm10")
ChIPQCreport(
    chip_obj,
    reportName = "ChIP_QC",
    reportFolder = "ChIPQC"
)
