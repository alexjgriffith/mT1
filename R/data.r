## This file is part of mT1
## http://github.com/alexjgriffith/mT1/, 
## and is Copyright (C) University of Ottawa, 2016. It is Licensed under 
## the GPL License; see LICENSE.txt.
## Author : Alexander Griffith
## Contact: griffitaj@gmail.com

#' IUPACDNA
#' 
#' A list with IUPAC names "A C G T R Y S W K M B D H V N" and
#' base pair values. eg IUPACDNA["A"] -> "A" &
#' IUPACDNA["N"]->c("A","C","T","G")
#' @format List with 15 members one for each IUPAC charcter :
"IUPACDNA"


#' Fasta
#' 
#' The nucleotide values for the genomic locations defined by peaks
#' @format A DANStringSet with 5670 members:
"mT1_fasta"

#' Jaspar
#'
#' A ist containing composite motifs and names for the jaspar database
#' @format A List with 4 members:
#' \describe{
#'   \item{jnames}{The composite motifs}
#'   \item{jfid}{The names of the moitfs, same order as jnames}
#'   \item{jl}{The numner of times a motif occured in the sample fasta file}
#'   \item{jsublM}{jnames sorted by jl}
#' }
"mT1_jaspar"

#' sampleMT1
#'
#' A sample mT1 object for use in examples.
#'
#' It was derived using mT1_fasta and the motifs
#' c("CANNTG","GATAA",mT1_jaspar$jsublM[1:10])
#' @format mT1 object:
#' 
"mT1_sampleMT1"

#' peaks
#'
#' A set of TAL1 peaks combined from 22 Hematopoietic and Endothelial
#' datasets
#'
#' @format A data frame with 5670 rows and 3 variables:
#' \describe{
#'   \item{chr}{chromosome}
#'   \item{start}{starting base pair location}
#'   \item{end}{ending base pair location}
#' }
"mT1_peaks"
