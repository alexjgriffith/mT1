# mT1
## Motifs in Tandem With One Another
### Alexander Griffith
### Aug 24 2016

**mT1** identifies the distance between DNA motifs in subsets of the genome. First it builds the empirical PDF from individual motifs, then it generates the expectations for the distances between the two motifs. **mT1** was designed to be applied to ChIP-Seq peaks.
___
### Installation
Currently only the source is available, compiled versions for Windows and Mac will be made available once mT1 is submitted to Bioconductor. There are two easy ways to install the package, via command line, and within R itself. The `devtools` package is required to install from within R.


From the command line:

```sh
## Clone the repository
git clone https://github.com/alexjgriffith/mT1.git
cd mT1

## Install 
R CMD INSTALL .
```
From within R:


```R
## Check if devtools is installed
if(!require(devtools)){
	install.packages("devtools")
	library(devtools)
}
## Install mT1
devtools::install_github("alexjgriffith/mT1")

## Install mT1 from the develop branch
## devtools::install_github("alexjgriffith/mT1",ref="develop")
```

The main branch is guaranteed to pass `R CMD check .` with no warnings once mT1 is build using `R CMD build --resave-data .`. The merges to the develop branch should always pass, however it is not a guarantee.


Note that mT1 also requires the Bioconductor package `Biostrings`.

To install `Biostrings`:

```R
source("http://bioconductor.org/biocLite.R")
biocLite("Biostrings")
```

___
### Examples
Compare the distances between two motifs in a set of genome subsets.

```R
## libraries
library(mT1)
## load the fasta file for analysis
fasta<-mT1_fasta

## Find the distances between CANNTG and GATAA for on each string
## they share. The distances are relative to the start of each
## motif
distances<-motifDistance(fasta,"CANNTG","GATAA")

## Plot the results between -20 and 20 bp
x<-seq(-20,20)
y<-combHeights(x,distances[,2])[[1]]
plot(x,y,main="CANNTG-GATAA",xlab="distance(bp)",
	 ylab="Frequency",type="l")

## Generate the emperical PDF
width<-300 # the width the fasta lines (must be uniform)
m1PDF<-motifPDF(fasta,"CANNTG")[,2]
m2PDF<-motifPDF(fasta,"GATAA")[,2]

## Find the expected probability for each distance
## Requires all peaks to be the same width
mp<-eMP(m1PDF,m2PDF,width)

## r is the range of the resulting mp, note that its length should 
## be one less than twice the width
r<-seq(-width+1,width-1)
hs<-combHeights(r,distances[,2])[[1]]
n<-length(distances[,2]) # total number of occurrences

## Apply the bionomial test to the expected probability and
## number of peaks seen
pvalue<-mapply(function(a,b)btest(a,n,b),hs,mp)

par(mfrow=c(2,1))
plot(r,pvalue)
plot(r,hs)

```


Analysis of multiple motifs under peaks that have been previously identified. This example requires `BSgenome.Hsapiens.UCSC.hg19`, which can be installed using `biocLite`.
```R
## libraries
library(mT1)
library(Biostrings) # needed to get fasta data from genomic co-ords
library(BSgenome.Hsapiens.UCSC.hg19) # genome

## load a set of Jaspar motifs as strings
jaspar<-mT1_jaspar

## Example set of peaks
peaks<-mT1_peaks

## Transform genomic co-ords into neucleotides
genome<-BSgenome.Hsapiens.UCSC.hg19
fasta<-getSeq(genome,peaks$chr,start=peaks$start,end=peaks$end)

## Motifs to compare
motifs<-c("CANNTG","GATAA",jaspar$jsublM[1:8])

## Find the preferred distances between `motifs` under the peaks
objMT1<-mT1(fasta,motifs)

## plot based on string
plot(objMT1,"CANNTG","GATAA")

## Plot based on printed order
print(objMT)
plot(objMT1,i=1)

```
