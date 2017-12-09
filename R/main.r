## This file is part of mT1
## http://github.com/alexjgriffith/mT1/, 
## and is Copyright (C) University of Ottawa, 2016. It is Licensed under 
## the GPL License; see LICENSE.txt.
## Author : Alexander Griffith
## Contact: griffitaj@gmail.com

#' Is IUPAC
#' 
#' Check if a string contains only iupac characters
#' @param string An atomic character
#' @return A logical, True if all IUPAC false otherwise
#' @examples
#' isIUPAC("CANNTG")
#' isIUPAC("CAXCTG")
#' @export
isIUPAC<-function(string){
    if(is.null(IUPACDNA))
        stop("missing IUPACDNA")
    spl<-strsplit(string,split="")[[1]]
    all(unlist(lapply(spl,function(s)!is.null(unlist(IUPACDNA[s])))))
}

#' IUPAC to Base
#'
#' This function Transforms an IUPAC DNA sequence into a string of 
#' neucleotides. Combo characters are replaced by their components in 
#' bracket. For example CNC -> C[ACGT]C. With the rl flag set to TRUE the
#' expanded set of sequences will be returned. For example CNC
#' would become c("CAC","CCC","CGC","CTC").
#' @param char the input string of IUPAC characters
#' @param rl flag to return all variants
#' @return Character containing neucleotides + brackets
#' @examples
#' mT1:::.IUPACtoBase("HGATAA")
#' mT1:::.IUPACtoBase("CANNTG")
#' mT1:::.IUPACtoBase("CANNTG",TRUE)#!/usr/bin/env R
.IUPACtoBase<-function (char, rl = FALSE){
    ## break atomic Character into a Character vector where the nchar
    ## in each Character is 1
    IUPAC <- strsplit(char, "")[[1]]
    if(is.null(IUPACDNA))
        stop("missing IUPACDNA")
    ## List of IUPAC characters and their composition
    ## IUPAC characters let you represent possible combinations of
    ## neucleotides as a single character. For example N = [ACGT]
    IUPACCharacters<-IUPACDNA
    ## Ensure all characters in char are IUPAC
    if (any(!IUPAC %in% names(IUPACCharacters))) {
        stop("Input string contains non IUPAC characters.")
    }
    ## Transform IUPAC characters to their component base pairs
    vals <- sapply(IUPAC, function(x) {
        (IUPACCharacters[x])
    })
    ## If the rl flag is set to true, set Base to a list of the
    ## component base pairs
    if (rl) {
        Base <- c("")
        for (i in vals) {
            kit <- c()
            for (j in Base) kit <- c(kit, paste(j, i, sep = ""))
            Base <- kit
        }
    }
    ## If the rl flag is false (default) concat all Base pair components.
    ## If the components are not singleton enclose them in brackets. This
    ## form allows for easy use in grep.
    else {
        Base <- c("")
        for (i in vals) if (length(i) == 1) 
            Base <- paste(Base, i, sep = "")
        else Base <- paste(Base, "[", do.call(paste, as.list(c(i, 
            sep = ""))), "]", sep = "")
    }
    return(Base)
}

#' Complement
#'
#' Takes a string of neucleotides (ACTG) and returns the complement
#' (TGAC). This works on strings only, for DNAStringSet please refer
#' to the `Biostrings::complement` function. To transform IUPAC characters
#' into neucleotide form refer to \code{\link{IUPACtoBase}}.
#' @param string the input string of nucleotides 
#' @return The compliment composite string as an atomic Character
#' @examples
#' mT1:::.complement(mT1:::.IUPACtoBase("CANNTG"))
#' mT1:::.complement("CA[AC]TT[ACGT]GG")
.complement<-function (string){
    ## Possible base pairs + Brackets
    chars <- c("A", "G", "C", "T", "[", "]")
    names(chars) <- c("T", "C", "G", "A", "]", "[")
    ## Generate the complement
    paste(rev(sapply(strsplit(string, "")[[1]], function(x) {
        (chars[x])
    })), collapse = "")
}

#' Get JASPAR
#'
#' A function to be used in conjunction with the mT1_jaspar.
#' Pass `getJASPAR` a composite motif from the jaspar database
#' and it will return the name of the motif.
#' @param name A composite motif. IUPAC DNA characters only.
#' @return The name of the motif as an atomic character, if it came
#' from the jaspar database, otherwise NULL.
#' @examples
#' # load the required Jaspar motif data
#' a<-mT1_jaspar$names[1]
#' getJASPAR(a)
#' @export
getJASPAR<-function(name){
    with(mT1_jaspar,{
        ## make sure required global variables are loaded
        if(is.null(mT1_jaspar$jfid)||is.null(mT1_jaspar$jnames)){
            stop(paste0("jaspar.RData must be loaded. ",
                        "load(system.file(\"extdata\",\"jaspar.RData\"",
                        ",package=\"mT1\"))"))
        }
        ## Find appropriate jasapar final id from the name index
        ret<-jfid[jnames == name]
        ## Test if null, currently this step is redundant
        ifelse(is.null(ret),NULL,ret)
    })
}

#' Motif Distance
#'
#' Find the distances between two motifs within set of genomic ranges.
#' The subsets are represented by DNAStringSet from the Biostrings package.
#' If the index of genomic ranges for each motif is already known use
#' \code{\link{diffMotif}}. If the distribution of motifs are already known
#' look into mT1:::.motifDiffFromDens.
#' 
#' @param fasta  A DNAStringSet with a set of genomic ranges
#' @param motif1 A motif as an atomic Character of IUPAC characters
#' @param motif2 A motif as an atomic Character of IUPAC characters
#' @return a numeric vector represented by two columns, the first being
#' the index the second being the position of the motif within that index.
#' distance
#' @examples
#' ## load fasta file
#' fasta<-mT1_fasta
#' ## Find distances
#' distances<-motifDistance(fasta,"CANNTG","HGATAA")
#'
#' ## Determine the range of the motif locations and plot the
#' ## frequency
#' width<-fasta@@ranges@@width[1] # all fasta strings should be the same width
#' r<-seq(-width+1,width-1)
#' y<-combHeights(r,distances[,2])[[1]]
#' plot(r,y,main="CANNTG-HGATAA")
#' @export
motifDistance<-function(fasta,motif1,motif2){
    motifs<-c(motif1,motif2)
    ## Find the location of all motifs
    mloc<-lapply(motifs,function(x) grep(.IUPACtoBase(x), fasta))
    ## Find the location of all motif complements
    cloc<-lapply(motifs,function(x) grep(.complement(.IUPACtoBase(x)),fasta))
    names(mloc)<-motifs
    names(cloc)<-motifs
    ## Call diff motifs to find the difference distribution between motif1
    ## and motif2
    diffMotif(fasta,motifs,mloc,cloc)
}

#' Motif Difference
#' 
#' Find the distances between two motifs within a subset of the genome
#' represented by fasta. If the distribution of motifs are already known
#' look into mT1:::.motifDiffFromDens.
#' @param fasta A DNAStringSet
#' @param motifs The list of motif Character Atomics
#' @param mloc A list of indicies  of each motif in fasta
#' @param cloc A list of indicies of each composite motif in fasta
#' @param i Motif index 1
#' @param j Motif index index 2
#' @param min numerical atomic, minimum length of mloc[[i or j]]
#' @param combine What to do if two motifs occur in the same region
#' Merged|First
#' @return A numeric vector with two columns , the first being the
#' genomic index the second being the distance.
diffMotif<-function(fasta,motifs,mloc,cloc,i=1,j=2,min=200,combine="Merged"){
    ## Select the motifs from index i an and j
    n1<-motifs[[i]]
    n2<-motifs[[j]]
    ## Find the genomic regions that are shared between motifs
    sharedm<-intersect(mloc[[n1]],mloc[[n2]])
    sharedc<-intersect(cloc[[n1]],cloc[[n2]])    
    both<-intersect(sharedc,sharedm)
    ## if there are not enough overlaping locations return NA
    if(length(union(sharedm,sharedc))<min)
        return(NA);
    ## Find the position of each motif in genomic regions that are shared
    l1m<-gregexpr(.IUPACtoBase(n1), fasta[sharedm,],ignore.case=TRUE)
    l2m<-gregexpr(.IUPACtoBase(n2), fasta[sharedm,],ignore.case=TRUE)
    l1c<-gregexpr(.complement(.IUPACtoBase(n1)),
                  fasta[sharedc,],ignore.case=TRUE)
    l2c<-gregexpr(.complement(.IUPACtoBase(n2)),
                  fasta[sharedc,],ignore.case=TRUE)
    ## Locations under the same peak are combined in two ways
    ## First: Keep the smallest distance
    ## Merged: Keep all distances   
    if(combine == "First"){
        ## Determine the distance between motifs
        mh<-cbind(loc=sharedm,pos=mapply(function(a,b){
            x<-unique(c(outer(a,b,Vectorize(function(i,j) i-j))))
            x[which.min(abs(x))]},l1m,l2m))
        ch<-cbind(loc=sharedc,pos=mapply(function(a,b){
            x<-unique(c(outer(a,b,Vectorize(function(i,j) j-i))))
            x[which.min(abs(x))]},l1c,l2c))
        ## Ensure that regions that have a motif and its complement found only
        ## return a single distance
        p2<-ch[!ch[,1] %in% both,]
        p3<-mh[!mh[,1] %in% both,]        
        ret<-rbind(p2,p3)
        if(length(both)>0){
            p1<-cbind(both,apply(cbind(ch[ch[,1] %in% both,2],
                                       mh[mh[,1] %in% both,2]),1,
                                 function(x) x[which.min(abs(x))]))
            ret<-rbind(ret,p1)
        }
    }
    else if(combine == "Merged"){
        ## A function that takes the distances between motifs and returns
        ## a two column numerical array with the distances and their
        ## genomic index.
        keep<-function(l1,l2,fun,shared){
            pos<-mapply(function(a,b) {
                unique(c(outer(a,b,Vectorize(fun))))},l1,l2,SIMPLIFY=FALSE)
            lens<-lapply(pos,length)
            loc<-unlist(mapply(function(x,n){
                rep(x,n)} ,shared,lens,SIMPLIFY=FALSE))
            cbind(loc=loc,pos=unlist(pos))
        }
        ## Difference in motif length;
        dchar<-nchar(n2)-nchar(n1)
        ## Determine the distances between motifs
        mh<-keep(l1m,l2m,function(i,j){i-j},sharedm)
        ch<-keep(l1c,l2c,function(i,j){j-i+dchar},sharedc)
        ## Make sure palandromic motifs are not listed twice
        p2<-ch[!ch[,1] %in% both,]
        p3<-mh[!mh[,1] %in% both,]        
        x<-merge(as.data.frame(ch),as.data.frame(mh))
        ret<-rbind(do.call(cbind,c(x[order(x[,1]),])),p2,p3)
    }        
    colnames(ret)<-c("loc","pos")
    ## A function used to make sure that occurances where the two motifs
    ## overlap are set to 0. This prevents the occurrence of high valued
    ## variants such as CANNTG and CACCTG with a prefered distance of 0.
    clean<-function(df){
        x<-df[,2]
        change<- x>= (- nchar(n2)) & x<= nchar(n1)        
        df[(!change), ]
    }
    clean(ret[order(ret[,1]),])
}

#' Motif PDF
#'
#' Generate a distribution of locations for a single motif within a set
#' of genomic regions.
#' @param fasta A DNAStringSet from Biostrings.
#' @param motif An atomic Character containing only DNA IUPAC characters
#' @return A numeric vector with two columns , the first being the index
#' the second being the distance
#' @examples
#' ## load fasta file
#' fasta<-mT1_fasta
#' pdf<-motifPDF(fasta,"CANNTG")
#' x<-seq(1:100)
#' y<-combHeights(x,pdf[,2])[[1]]
#' ## plot the results
#' plot(x,y,main="CANNTG",ylab="Frequency",xlab="Index")
#' @export
motifPDF<-function(fasta,motif){
    ## Find the indicies of the genomic ranges where the motif and
    ## motif complement are present
    mloc<-grep(.IUPACtoBase(motif), fasta)
    cloc<-grep(.complement(.IUPACtoBase(motif)),fasta)
    ## Use findlocs to build the distribution.
    findLocs(fasta,mloc,cloc,motif)
}

#' Find Motif Locations
#'
#' Finds the position of motifs on each index of fasta. The results can
#' be used to build the emperical PDF of motif locations.
#' @param fasta A DNAStringSet
#' @param mloc fasta The indicies that have  motifs
#' @param cloc fasta The indicies that have compliment motifs
#' @param n1 The motif name, an Atomic Character of only IUPAC characters
#' @param combine If two motifs are found on the same index Merged|First
#' @return Two columns , the first being the instance the second being the
#' distance
findLocs<-function(fasta,mloc,cloc,n1,combine="Merged"){
    ## Find the positions of the individual motif n1
    l1m<-gregexpr(.IUPACtoBase(n1), fasta[mloc],ignore.case=TRUE)
    l1c<-gregexpr(.complement(.IUPACtoBase(n1)),fasta[cloc],ignore.case=TRUE)
    both<-intersect(mloc,cloc)
    ## The remainder of this function is very similar to motifDiffs
    ## applied to a single motif
    keep<-function(l,loc){
        pos=unlist(lapply(l,unique))
        len<-sapply(l,length)
        eloc<-unlist(mapply(function(x,n)rep(x,n),loc,len,SIMPLIFY=FALSE))
        cbind(loc=eloc,pos=pos)
    }
    if(combine == "First"){
        mh<-cbind(loc=mloc,sapply(l1m,function(x) x[which.min(abs(x))]))
        ch<-cbind(loc=cloc,sapply(l1c,function(x) x[which.min(abs(x))]))    
        p2<-ch[!ch[,1] %in% both,]
        p3<-mh[!mh[,1] %in% both,]
        ret<-rbind(p2,p3)
        if(length(both)>0){        
            p1<-cbind(both,apply(cbind(ch[ch[,1] %in% both,2],
                                       mh[mh[,1] %in% both,2]),1,
                                 function(x) x[which.min(abs(x))]))
            ret<-rbind(ret,p1)
        }
        colnames(ret)<-c("loc","pos")
    }
    else if(combine == "Merged"){        
        mh<-keep(l1m,mloc)
        ch<-keep(l1c,cloc)
        p2<-ch[!ch[,1] %in% both,]
        p3<-mh[!mh[,1] %in% both,]
        x<-merge(as.data.frame(ch),as.data.frame(mh))
        ret<-rbind(do.call(cbind,c(x[order(x[,1]),])),p2,p3)
        colnames(ret)<-c("loc","pos")
    }
    else{
        ## This option should be used if you want to sort through the
        ## occurances of motifs and their complements later.
        mh<-cbind(keep(l1m,mloc),0)
        ch<-cbind(keep(l1c,cloc),1)
        ret<-rbind(mh,ch)
        colnames(ret)<-c("loc","pos","dir")
    }        
    ret[order(ret[,1]),]        
}

#' Make Title
#'
#' A generic function that  Generates a title to be used for plotting
#' an object.
#' @param x The obj to make title from
#' @param ... possible properties
#' @rdname makeTitle
makeTitle<-function(x,...){
    UseMethod("makeTitle",x)
}


#' @return An Atomic Character with motif info
#' @param i The index of interest from print(x)
#' @rdname makeTitle
#' @method makeTitle mT1
#' @export
makeTitle.mT1<-function(x,i,...){
    with(x,{
        do.call(paste,append(lapply(combs[i,],function(x)
            motifs[x]),list(sep="-")))
    })
}

#' @method summary mT1
#' @export
summary.mT1<-function(object,...){
    main<-data.frame(t(apply(object$combs[!object$sig,],1,
                             function(i) object$motifs[i])))
    xa<-seq(-object$width+1,object$width-1)
    if(length(main)>0){
        colnames(main)<-c("motif1","motif2")
        mpv<-sapply(which(!object$sig),function(i)
            min(object$pvalue[[i]]))
        mpvl<-sapply(which(!object$sig),function(i)
            xa[which.min(object$pvalue[[i]])])
        return(data.frame(main,mpv,mpvl))
    }
    else
        return (NULL)
}

#' @method print mT1
#' @export
print.mT1<-function(x,...){
    #cat("mT1\n\n")
    cat("motifs:  ")
    cat(paste0(x$motifs[1],"\n"))
    cat(paste("        ",x$motifs[2:length(x$motifs)],collapse="\n"))
    cat("\n\n")
    cat(paste0("combs: ",length(x$combs), "\n"))
    cat(paste0("sufficent: ",sum(!x$sig)),"\n")
    main<-summary.mT1(x)
    print(main)
}

#' @method plot mT1
#' @export
plot.mT1<-function(x,motif1=NULL,motif2=NULL,i=NULL,...){
    if(!is.null(i)){
        main<-data.frame(t(apply(x$combs[!x$sig,],1,function(i) x$motifs[i])))
        motif1<-as.character(main[i,1])
        motif2<-as.character(main[i,2])
    }
    with(x,{        
        xa<-seq(-width+1,width-1)
        n1<-which(motifs == motif1)[1]
        n2<-which(motifs == motif2)[1]
        i<-union(which(combs[,1] == n1 & combs[,2] == n2),
                 which(combs[,1] == n2 & combs[,2] == n1))[1]
        par(mfrow=c(3,2))
        plot(density(dens[[n1]][,2]-width/2),main=motif1)
        plot(density(dens[[n2]][,2]-width/2),main=motif2)
        plot(xa,mp[[i]],type="l",ylab="p",xlab="Index",main="Convolution")
        plot(xa,hs[[i]],type="l",xlab="Index",ylab="Frequency",
             main=paste0(motifs[combs[i,1]],"-",motifs[combs[i,2]]),...)
        plot(combHeights(seq(0,max(hs[[i]])),
                         hs[[i]])[[1]],type="l",xlab="Height",
             ylab="Frequency",main="Freq")
        plot(xa,pvalue[[i]],type="l",ylab="p-value",xlab="Index",
             main="Tests",...)
    })
}

#' Motif Difference From Density
#'
#' If two set of motif locations and positions are known this function
#' can be used to determine the difference in motif location.
#' 
#' @param a The distribution of motif a
#' @param b The distribution of motif b
#' @param diff The difference in motif size between a and b
#' @return A two column vector with location indexes and positions
.motifDiffFromDens<-function(a,b,diff){
    ## todo: clean up .motifDiffFromDens, giving better names to t2 and t3
    af<-as.data.frame(a)
    bf<-as.data.frame(b)
    colnames(af)<-c("loc","a","dir")
    colnames(bf)<-c("loc","b","dir")
    t2<-merge(af[af$dir == 1,c(1,2)],bf[bf$dir == 1,c(1,2)])
    back<-data.frame(loc=t2[,1],pos=t2[,3]-t2[,2]+diff)
    t3<-merge(af[af$dir == 0,c(1,2)],bf[bf$dir == 0,c(1,2)])
    forw<-data.frame(loc=t3[,1],pos=t3[,2]-t3[,3])    
    ret<-rbind(back,forw)    
    #ret<-rbind(cbind(t3[,1],forw,0),cbind(t2[,1],back,1))
    #colnames(ret)<-c("loc","pos","dir")
    ret[!duplicated(ret),]
}



#' Motifs in Tandem with One Another
#'
#' Generates a mT1 object from DNAStringSet and a vector of motifs. A mT1
#' object can be used for the investigation of preferred motif positions
#' and identifying composite motifs. There must be at least 3 motifs being
#' compared. For the comparison of two motifs look at motifDiff. The purpose
#' of a mT1 object is to identify novel preferred distances from a large
#' set of motifs.
#' @param fasta The set of genomic locations as a DNAStringSet
#' @param motifs The list of three or more motifs being compared,
#' IUPAC chars only.
#' @param verbose A flag to print motifs as completed.
#' @param cl A cluster from makeForkCluster
#' @return A mT1 object
#' @examples
#' library(Biostrings) # needed to get fasta data from genomic co-ords
#' library(BSgenome.Hsapiens.UCSC.hg19) # genome
#' ## load a set of Jaspar motifs as strings
#' ## load(system.file("extdata","jaspar.RData",package="mT1"))
#' jaspar<-mT1_jaspar
#' ## Example set of peaks
#' ## load(system.file("extdata","peaks.RData",package="mT1"))
#' peaks<-mT1_peaks
#' ## Transform genomic co-ords into neucleotides
#' genome<-BSgenome.Hsapiens.UCSC.hg19
#' fasta<-getSeq(genome,peaks$chr,start=peaks$start,end=peaks$end)
#' 
#' ## Motifs to compare
#' motifs<-c("CANNTG","GATAA",jaspar$jsublM[1:2])
#' 
#' ## Find the preferred distances between `motifs` under the peaks
#' my_sampleMT1<-mT1(fasta,motifs)
#' 
#' ## plot based on string
#' plot(my_sampleMT1,"CANNTG","GATAA")
#' 
#' ## Plot based on printed order
#' print(my_sampleMT1)
#' plot(my_sampleMT1,i=1)
#' @export
mT1<-function(fasta,motifs,verbose=FALSE,cl=NULL){
    makePind<-function(l,n,name=NULL){
        if(!is.null(name))
            cat(name)
        if(l<n){
            pind<-seq(l)
            per <- 1/(l-2)*100
            pb <- txtProgressBar(style=3,width=60,title=name)
        }
        else{
            pind<-c(ceiling(seq(1,l-ceiling(l/(n+1)),by=l/(n+1))),l)
            per<-1/n*100
            pb <- txtProgressBar(style=3,width=60,label=name)
        }
        return(list(pind=pind,per=per,pb=pb,go=TRUE,name=name))
    }
    checkPind<-function(complete,pind){
        if(pind$go){
            perc<-(which(complete ==pind$pind)-1)*pind$per
            if(length(perc)==1){                    
                setTxtProgressBar(pind$pb, perc/100,label=pind$name)
                if(which(complete==pind$pind)==length(pind$pind) |
                   perc/100==1){
                    close(pind$pb)
                    pind$go<-FALSE
                }

            }
        }
        pind
    }
    ## Make sure motifs are unique
    tmots<-unique(motifs)
    if(length(tmots)!=length(motifs)){
        warning("motifs must be unique")
        motifs<-tmots
    }
    ## Ensure that there are 3 motifs at least. For the comparison of
    ## two motifs refer to motifDiff    
    if(length(motifs)<3){
        warning("must provide mT1 with more than 2 motifs")
        return (NA)
    }
    ## find which fasta indicies have the motifs of interest
    if(verbose){
        pind<-makePind(length(motifs),20,"Finding motif co-ords ...\n")
    }
    mloc<-lapply(motifs,function(x){
        if(verbose){
            pind<-checkPind(which(motifs==x),pind)
        }
        grep(.IUPACtoBase(x), fasta)})
    if(verbose){
        pind<-makePind(length(motifs),20,"Finding complement co-ords ...\n")
    }
    cloc<-lapply(motifs,function(x){
        if(verbose){
            pind<-checkPind(which(motifs==x),pind)
        }
        grep(.complement(.IUPACtoBase(x)),fasta)})
    names(mloc)<-motifs
    names(cloc)<-motifs
    ## Determine the individual PDFs for each motif
    if(verbose){
        pind<-makePind(length(motifs),20,"Finding motif locations ...\n")
    }
    a<-lapply(motifs,function(x){
        if(verbose){
            pind<-checkPind(which(motifs==x),pind)
        }
        findLocs(fasta,mloc[[x]],cloc[[x]],x,"FALSE")})
    ## combination of all motifs
    tofind<-t(combn(seq(length(motifs)),2))

    if(verbose){
        pind<-makePind(dim(tofind)[1],20,
                       paste0("Difference Analysis Begining.",
                              " Searching through ",dim(tofind)[1],
                              " combinations\n"))        
    }
    t1<-apply(tofind,1,function(x){        
        if(verbose){
            checkPind(which(tofind[,1]==x[1] &  tofind[,2]==x[2]),pind)
        }
        m1<-a[[x[1]]]
        m2<-a[[x[2]]]
        diff<-nchar(motifs[x[2]])-nchar(motifs[x[1]])
        if(length(intersect(m1[,1],m2[,1]))>300){
            y<-.motifDiffFromDens(m1,m2,diff)        
            change<- y[,2]>= (- nchar(motifs[x[2]])) &
                y[,2]<= nchar(motifs[x[1]])
            y<-y[(!change), ]
        }
        else
            y<-NA
        y
    })
    ## Which indicies have values
    large<-sapply(t1,function(x) all(unlist(is.na(x))))
    mT1<-list(diff=t1,dens=a,sig=large,motifs=motifs,combs=tofind)
    ## determine the p-values for each motif comb
    ## If biostrings is loaded nchar(fasta) will be equivalent
    ## Check this after any update to IRanges
    width<-fasta@ranges@width[1]
    ## Determine probabilities
    if(verbose){
        pind<-makePind(dim(tofind)[1],20, "Finding P-Values ...\n")
    }
    prob<-apply(cbind(tofind,seq(dim(tofind)[1])),1,
                function(x){
                    if(verbose)
                        checkPind(x[3],pind)                    
                    .ePD(t1[[x[3]]],a[[x[1]]][,2],a[[x[2]]][,2],
                                 width,motifs[[x[2]]])})
    if(verbose)
        cat("Finalizing mT1 object ...\n")
    ## Build and return the mT1 object
    mT1<-append(mT1,
                list(fasta=fasta,width=width,hs=lapply(prob,function(x) x[["hs"]]),
                mp=lapply(prob,function(x) x[["mp"]]),
                pvalue=lapply(prob,function(x) x[["pvalue"]])))
    ## Set class type mT1
    attr(mT1,"class")<-"mT1"
    mT1
}

#' @method c mT1
#' @export
c.mT1<-function(...,recursive=FALSE){
    wl<-list(...)
    if(length(wl)<2)
        return (wl[[1]])
    ## this function can be greatly simplified by looping over
    ## c rbind and first
    toRbind<-c("combs")
    toC<-c("diff","dens","sig","motifs","hs" ,"mp","pvalue")
    toTest<-c("fasta","width")
    res<-list()
    for(n in toC){        
        res[[n]]<-do.call(c,lapply(wl,function(x){
            if(is.null(x[[n]]))
                stop(paste0(n," not found in mT1"))
            x[[n]]}))
    }
    for(n in toRbind){
        res[[n]]<-do.call(rbind,lapply(wl,function(x){
            if(is.null(x[[n]]))
                stop(paste0(n," not found in mT1"))
            x[[n]]}))
    }
    fasta<-lapply(wl,function(x) x[["fasta"]])
    if(!Reduce(function(a,b)a == all(b == fasta[[1]]),
               init=TRUE,x=fasta))
        stop("different genomic ranges")
    width<-unlist(lapply(wl,function(x) x[["width"]]))
    if(any(!(width[1] == width)))
        stop("different widths")
    else{        
        width <- width[1]
        
    }
    first<-function(...){
        lapply(list(...),function(x) x[[1]])
    }
    for(n in toTest){
        res[[n]]<-do.call(first,lapply(wl,function(x){
            if(is.null(x[[n]]))
                stop(paste0(n," not found in mT1"))
            x[[n]]}))
    }
    res$width<-width
    attr(res,"class")<-"mT1"
    return( res)
}

#' Add Motif
#'
#' Takes a motif and creates a new object based on the class of object
#'
#' @param object A object with an addMotif method
#' @param motif A Atomic Character containing IUPAC characters
#' @param ... Variables to pass to object.
#' @return mT1 object with new motif
#' @rdname addMotif
#' @examples
#' new<-addMotif(mT1_sampleMT1,"CACCTG")
#' @export
addMotif<-function(object,motif,...){
    UseMethod("addMotif",object)
}

#' @rdname addMotif
#' @method addMotif default
#' @export
addMotif.default<-function(object,motif,...){
    motif
}

#' @rdname addMotif
#' @method addMotif mT1
#' @export
addMotif.mT1<-function(object,motif,...){
    ## Much of this function is shared with mT1
    ## I need to extract the similarities and move them
    ## to a set of private mT1 helper functions
    return ( with(object,{
        mloc<-grep(.IUPACtoBase(motif), fasta)
        cloc<- grep(.complement(.IUPACtoBase(motif)),fasta)        
        a<-findLocs(fasta,mloc,cloc,motif,"FALSE")
        nm1<-length(motifs)
        tofind<-cbind(nm1+1,seq(nm1))
        t1<-apply(tofind,1,function(x){
                m1<-a
                m2<-dens[[x[2]]]
                diff<-nchar(motifs[x[2]])-nchar(motif)
                if(length(intersect(m1[,1],m2[,1]))>300){
                    y<-.motifDiffFromDens(m1,m2,diff)        
                    change<- y[,2]>= (- nchar(motifs[x[2]])) &
                        y[,2]<= nchar(motif)
                    y<-y[(!change), ]
                }
                else
                    y<-NA
                y
        })
        prob<-apply(cbind(tofind,seq(dim(tofind)[1])),1,
                    function(x) .ePD(t1[[x[3]]],a[,2],dens[[x[2]]][,2],
                                     width,motifs[[x[2]]]))
        large = sapply(t1,function(x) all(unlist(is.na(x))))
        ret<-list(diff=t1,dens=list(a),sig=large,motifs=motif,combs=tofind,
             width=width,
             fasta=fasta,
             hs=lapply(prob,function(x) x[["hs"]]),
             mp=lapply(prob,function(x) x[["mp"]]),
             pvalue=lapply(prob,function(x) x[["pvalue"]]))
        class(ret)<-"mT1"
        c(object,ret)
    }))
}

#' Binomial Test
#'
#' A wrapper around binom.test that ensures that n!=0 and that k is at
#' least a minimum number.
#' @param k count
#' @param n total
#' @param p prob
#' @param kmin minimum k cut off
#' @param nmin minimum n cut off
#' @return a atomic numeric p-value
#' @examples
#' ## How many times are two motifs found D base pairs apart
#' frequencyAtD<-10
#' ## How may times are two motifs found any D appart
#' possibleAtD<-200
#' ## If one motif distance was computed what is the probability of it
#' ## being D
#' probabilityAtD<-0.03
#' btest(frequencyAtD,possibleAtD,probabilityAtD)
#' @export
btest<-function(k,n,p,kmin=10,nmin=600){
    if(n<nmin) # n must be greater than 0
        return(0)
    ## binom.test needs to be speed up, half of the time is spent here
    else if(k>kmin){ # ensure that there are more k than min
        ##return(log(binom.test(k,n,p)$p.value,10))
        return (log(dbinom(k,n,p),10))}
    else
        return(0)
}

#' Expected Motif Probability
#'
#' Determine the expected probability of finding two motifs at a specific
#' distance based on the empirical PDFs of the individual motifs
#' @param a numerical vector PDF of motif a
#' @param b numerical vector PDF of motif b
#' @param width width of fasta indicies
#' @param nb Difference in length between motif a and motif b
#' @return expectation for each motif distance 
#' @export
eMP<-function(a,b,width,nb=0){
    refl<-function(x,width,sizex){        
        abs(x-width+sizex-2)
    }
    y<-combHeights(seq(1,width),a,refl(b,width,nb))
    mod<-convolve(y[[1]],y[[2]],type="open")
    mod[mod<0]<-0
    mod/sum(mod)
}

#' Expected Probability Distance
#'
#' Determine the pvalue from a vector of distances and PDFs from the
#' two motifs being compared.
#' @param t1 A width 2 numerical vector of indicies and distances
#' @param a A numerical vector PDF of motif a
#' @param b A numerical vector PDF of motif b
#' @param width A numeric atomic, the width of fasta indicies, should
#' @param ma String motif A
#' @param mb String motif B
#' be uniform
#' @return list(hs,mp,p-values)
.ePD<-function(t1,a,b,width,mb){
    ## make sure there are rows in t1
    if(any(is.na(t1))){        
        return(list(hs=NA,mp=NA,pvalue=NA))
    }
    hs<-combHeights(seq((-1* width+1),width-1),t1[,2])[[1]]
    n<-sum(hs)
    if(max(hs)>n){
        stop("n >hs")
    }
    mp<-eMP(a,b,width,nchar(mb))
    pvalue<-dbinom(hs,n,mp,log=TRUE)
    pvalue[hs<1]<-0
    list(hs=hs,mp=mp,pvalue=pvalue)
}

#' Combination Heights
#'
#' Determines the number of occurrences in each member of ... along
#' sequence x.
#' @param x A sequence of possible values represented as indicies
#' @param ... The vectors of values to be mapped to x and summed
#' @return List with length of x range
#' @examples
#' x<-seq(20)
#' y<-runif(200,1,20)
#' par(mfrow=c(1,2))
#' plot(y)
#' plot(x,combHeights(x,y)[[1]])
#' @export
combHeights<-function(x,...){
    if(length(as.list(...)) == 0){
        warning("combHeights should be passed one additional variable")
        return(NA)
    }
    min<-min(x)
    max<-max(x)
    lapply(list(...),function(h){        
        y<- rep(0,length(x))
        for(i in h){
            if(i>=min & i<=max)
                y[i-min +1]<-y[i-min +1]+1
        }
        y
    })
}
