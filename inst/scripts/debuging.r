sample<-mT1_sampleMT1



names(sample)

a<-sample$dens[[1]][,2]

b<-sample$dens[[2]][,1]

eMP(a,b,300)

plot(sample,i=2)

ePD(sample$diff[[1]],sample$dens[[1]],sample$dens[[2]],300)



    
with(sample,{
    combHeights(seq(0,max(hs[[1]])),
                hs[[1]])
    })


mT1_fasta

sample<-mT1(mT1_fasta,c("CANNTG","GATAA",mT1_jaspar$jsublM[1:10]))



library(mT1)

                                        #sample<-mT1_sampleMT1

sample3<-mT1(mT1_fasta,c("CANNTG","GATAA",mT1_jaspar$jsublM[1:2]))

a<-sample3$dens[[1]]
b<-sample3$dens[[2]]


##a[,1] b[,1]


.motifDiffFromDens<-function(a,b,diff){
    af<-as.data.frame(a)
    bf<-as.data.frame(b)
    colnames(af)<-c("loc","a","dir")
    colnames(bf)<-c("loc","b","dir")
    t2<-merge(af[af$dir==1,c(1,2)],bf[bf$dir==1,c(1,2)])
    back<-data.frame(loc=t2[,1],pos=t2[,3]-t2[,2]+diff)
    t3<-merge(af[af$dir==0,c(1,2)],bf[bf$dir==0,c(1,2)])
    forw<-data.frame(loc=t3[,1],pos=t3[,2]-t3[,3])    
    ret<-rbind(back,forw)    
    #ret<-rbind(cbind(t3[,1],forw,0),cbind(t2[,1],back,1))
    #colnames(ret)<-c("loc","pos","dir")
    ret[!duplicated(ret),]
}


plot(combHeights(seq(-30,30),motifDiff(a,b,-1)[,2])[[1]])


library(mT1)

obj<-mT1(mT1_fasta,c("CANNTG","GATAA","HGATAA"))

plot(obj,i=2)


clean<-function(x,n1,n2){
        change<- x>= (- nchar(n2)) & x<= nchar(n1)        
        x[(!change) ]
    }


x<-seq(-15,15)
y<-combHeights(x,clean(.motifDiffFromDens(a,b,-1)[,2],"CANNTG","GATAA"))[[1]]
plot(x,y)




Rprof("~/Desktop/rprof.out")
sample3<-mT1(mT1_fasta,c("CANNTG","GATAA",mT1_jaspar$jsublM[1:2]))
summaryRprof("~/Desktop/rprof.out")

log(binom.test(250,12000,0.003)$p.value,2.7)

log(dbinom(250,12000,0.003),10)

sample3<-mT1(mT1_fasta,unique(c("CANNTG","GATAA",mT1_jaspar$jsublM[1:2])))

new<-addMotif(sample3,"CANNTG")

c(sample3,new)


library(CCCA)
library('BSgenome.Hsapiens.UCSC.hg19')
source("~/r-workspace/project/ccca.r")
source("~/r-workspace/project/project.r")
source("~/r-workspace/project/project-variables.r")

env<-getPRC20(2)
env<-addFasta(env,width=150)

devtools::install_github("alexjgriffith/mT1",ref="develop")

library("mT1")

largeAnalysis<-mT1(env$fasta,unique(c("CANNTG","GATAA",
                                      mT1_jaspar$jsublM[1:300])),verbose=TRUE)


largeAnalysis<-mT1(mT1_fasta,unique(c("CANNTG","GATAA",
                                      mT1_jaspar$jsublM[1:10])),verbose=TRUE)

cbind(largeAnalysis$hs[[1]],largeAnalysis$mp[[1]])

par(mfrow=c(1,2))
plot(combHeights(seq(1,150),abs(largeAnalysis$dens[[1]][,2]-151))[[1]])
plot(combHeights(seq(1,150),largeAnalysis$dens[[1]][,2])[[1]])


testit <- function(x = sort(runif(20)), ...)
{
    pb <- txtProgressBar(...)
    for(i in c(0, x, 1)) {Sys.sleep(0.5); setTxtProgressBar(pb, i)}
    Sys.sleep(1)
    close(pb)
}

testit(style = 3,width=60)


refl<-function(v,width,center=width/2){
    ##abs(v-width-1)
    x<-seq(width)
    <-cbind(v,width)
}


ebox<-combHeights(seq(1,150),c(floor(rnorm(10000,75,30)),floor(runif(10000,1,152))),c(floor(runif(10000,1,152))))


a<-ebox[[1]]/sum(ebox[[1]])
b<-ebox[[2]]/sum(ebox[[2]])
d<-convolve(b,c,type="open")

par(mfrow=c(2,2))
plot(a,type="l")
plot(b,type="l")
plot(rev(a),type="l")
plot(d,type="l")



b<-c(floor(runif(10000,1,146)))

par(mfrow=c(2,1))
plot(combHeights(seq(150),b)[[1]],type="l")
plot(combHeights(seq(150),abs(b-150-1+5))[[1]],type="l")


     
  
log.env$btest<-list()

lapply(log.env$btest,function(x) dbinom(x[,1],x[,2],x[,3]))



do.call(convolve,append(combHeights(seq(1,150),largeAnalysis$dens[[1]][,2], abs(largeAnalysis$dens[[2]][,2]-151)),list(type="open")))


a<-summary(largeAnalysis)

a[order(a[,3])[1:100],]

png("~/Desktop/GATA-HLTF.png")
plot(largeAnalysis,i=2)
dev.off()


png("~/Desktop/sox10-HMBOX1.png")
plot(largeAnalysis,i=1721)
dev.off()

238              GATAA           CATATGK  -128.30633   15


2517          RTCTGGHW           TTTTCCA  -129.56547   13

png("~/Desktop/Tcf3-NFATC2.png")
plot(largeAnalysis,i=2517)
dev.off()


png("~/Desktop/Ebox-Gata.png")
plot(largeAnalysis,i=1)
dev.off()
