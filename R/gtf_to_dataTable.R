#' Extract attributes from GTF file
#'
#' @param gtf_attributes single line of GTF attribute column, e.g. genes$attributes[1]
#' @param att_of_interest string indicating which attribute value should be extracted, one of: c("gene_id", "transcript_id","")
#' @details Splits the attribute field of a GTF file line and returns just the value for the
#' attribute of interest.
#' Example usage: 
#' genes <- fread("~/Documents/Projects/2016-08_Walsh/data/gencode.v23.basic.lincRNA.gtf")
#' names(genes) <- c("chr","source","type","start","end","score","strand","phase","attributes")
#' genes$score <- NULL
#' genes$phase <- NULL
#' setkey(genes, chr, start, end)
#' pir.peaks.genes <- foverlaps(pir.peaks.merged, genes, nomatch = 0)
#' pir.peaks.genes$gene_id <- lapply(pir.peaks.genes$attributes, extract_attributes, "gene_id") %>% unlist
extract_attributes <- function(gtf_attributes, att_of_interest){
	library(data.table)
	att <- strsplit(gtf_attributes, "; ")
	att <- gsub("\"","",unlist(att))
	if(!is.null(unlist(strsplit(att[grep(att_of_interest, att)], " ")))){
    return( unlist(strsplit(att[grep(att_of_interest, att)], " "))[2])
  }else{
    return(NA)}
}
