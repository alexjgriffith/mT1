sample<-mT1_sampleMT1



names(sample)

a<-sample$dens[[1]][,2]

b<-sample$dens[[2]][,1]

eMP(a,b,300)

plot(sample,i=2)

ePD(sample$diff[[1]],sample$dens[[1]],sample$dens[[2]],300)


sample<-
    
with(sample,{
    combHeights(seq(0,max(hs[[1]])),
                hs[[1]])
    })


mT1_fasta

sample<-mT1(mT1_fasta,c("CANNTG","GATAA",mT1_jaspar$jsublM[1:10]))
