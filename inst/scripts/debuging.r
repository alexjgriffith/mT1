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

a<-sample2$dens[[1]]
b<-sample2$dens[[2]]


##a[,1] b[,1]


motifDiff<-function(a,b,diff){
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













