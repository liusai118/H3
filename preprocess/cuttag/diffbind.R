library(DiffBind)
abc <- read.csv("dataset.csv", stringsAsFactors = FALSE)
H3 <- dba(sampleSheet = abc)
mask4  <- dba.count(H3, minOverlap = 2,summits = 250)
info <- dba.show(mask4)
libsizes <- cbind(
  LibReads = info$Reads,
  FRiP = info$FRiP,
  PeakReads = round(info$Reads * info$FRiP)
)
rownames(libsizes) <- info$ID
print(libsizes)
peaks4 <- dba.peakset(mask4, bRetrieve = TRUE)
length(peaks4)
library(rtracklayer)
export(peaks4,
       con = "H3_minOv2_summits250_consensusPeaks.bed",
       format = "BED")

