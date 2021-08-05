# Section 2.2 : custom script to only keep read pairs where 
#one of the reads has two or more T nucleotides at the 5' end

#######################################################
##  Read me      ######################################
## Poly T criteria
## This function uses the package ShortRead to
## grab all the reads that start with two or more
## T's at the 5' end of either read in a pair.
## This corresponds to reads that are from the
## poly T tail of the mRNA 
##
##  Read me      ######################################
#######################################################

#######################################################
##  Dependencies      #################################

library(ShortRead)

##  Dependencies      #################################
#######################################################

setwd('/home/jhhal/regulation_de_novo/data/RNAseq/raw/schraiber/')

#######################################################
##  Filters and trims     #############################

fnFilterTT <- srFilter(function(x) {
  ## This filters out all reads that do not either
  ## start with at least two T's or at least two A's.
  ## 
  grepl('^TT', gsub('N', '', sread(x)))
},
name='FilterTT'
)

fnRunFilter <- function(read.files) {
  ## read.files is a vector with two file names,
  ## one for the fastq file with read 1 and
  ## another with the fastq file with read 2
  
  stopifnot(length(read.files) == 2)
  
  ## Read in the fastq files,
  ## reads is a list with two
  ## Fastq files, one for each
  ## read in a pair 
  reads <- lapply(read.files, readFastq)
  
  ## Run AATT filter
  ve.filterTT <- lapply(reads, fnFilterTT)
  
  ## Combine the the results from the filters
  ## on the reads of the different pairs in order
  ## to keep the reads that satisfy the condition
  ## in either of the pairs.
  ve.filterTT.paired <- ve.filterTT[[1]] | ve.filterTT[[2]]
  
  reads <- lapply(reads, function(x) x[ve.filterTT.paired])
  
  ## Write the files
  
  writeFastq(object = reads[[1]], compress = TRUE, mode = 'w',
             file = gsub('.fastq.gz', '_polyT-filtered.fastq.gz',
                         read.files[1]))
  
  writeFastq(object = reads[[2]], compress = TRUE, mode = 'w',
             file = gsub('.fastq.gz', '_polyT-filtered.fastq.gz',
                         read.files[2]))
  
  ## gc() removes the memory allocation
  gc()
}

##  Filters and trims     #############################
#######################################################

#######################################################
##  Running filter for all files   ###########

ve.files <- list.files(pattern = 'fastq')

li.files <- split(x = ve.files, f = gsub('-.*', '', ve.files))

lapply(li.files, fnRunFilter)

##  Running filter for all files   ###########
#######################################################
