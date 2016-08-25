library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)

## From the CCCA
load("~/Dropbox/UTX-Alex/Paper/Raw Data/peakLocations.RData")

mT1_fasta<-PCA$fasta
mT1_peaks<-PCA$bed
mT1_jaspar<-jaspar
mT1_sampleMT1<-mT1(mT1_fasta,c("CANNTG","HGATAA",mT1_jaspar$jsublM[1:10]))

save(mT1_fasta,file="mT1_fasta.RData")
save(mT1_peaks,file="mT1_peaks.RData")
save(mT1_jaspar,file="mT1_jaspar.RData")
save(mT1_sampleMT1,file="mT1_sampleMT1.RData")

