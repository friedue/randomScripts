#' @title Turn GRanges into data.table
#' @peaks.granges subset of GRanges object, preferably the result of RIPSeeker (RIPGRList)
#' @return data.table with the genomic coordinates and statistical measures per peak
ranges2dt <- function(peaks.granges){
  # As with IRanges objects, the accessors 
  #`start`, `end`, and `width` work, 
  # as well as some new ones: `seqnames`, `strand`, and `mcols`. 
  #`mcols` returns the DataFrame of metadata.

 library(data.table)

  a <- as.data.table(seqnames(peaks.granges)) 
  names(a) <- "chr"
  b <- as.data.table(ranges(peaks.granges))
  c <- as.data.table(mcols(peaks.granges))
  return(cbind(a,b,c))
}

