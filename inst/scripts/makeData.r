library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)

## From the CCCA
load("~/Dropbox/UTX-Alex/Paper/Raw Data/peakLocations.RData")

mT1_fasta<-PCA$fasta
mT1_peaks<-PCA$bed

rownames(mT1)<-NULL
    
mT1_jaspar<-jaspar


mT1_sampleMT1<-mT1(mT1_fasta,unique(c("CANNTG","GATAA",mT1_jaspar$jsublM[1:2])))

save(mT1_fasta,file="mT1_fasta.RData")
save(mT1_peaks,file="mT1_peaks.RData")
save(mT1_jaspar,file="mT1_jaspar.RData")

save(mT1_sampleMT1,file="mT1_sampleMT1.RData")

IUPACDNA <- list("A", "C", "G", "T", c("A", "G"), 
                        c("C", "T"), c("C", "G"), c("A", "T"),
                        c("G", "T"), c("A", 
                                                                              "C"), c("C", "G", "T"), c("A", "G", "T"), c("A", 
            "C", "T"), c("A", "C", "G"), c("A", "C", "G", "T"))
names(IUPACDNA) <- c("A", "C", "G", "T", "R", "Y", 
                            "S", "W", "K", "M", "B", "D",
                            "H", "V", "N")


save(IUPACDNA,file="IUPACDNA.RData")
