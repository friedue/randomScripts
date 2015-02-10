### run these lines once you start the R session
options(stringsAsFactors=FALSE)

# define your input data
FN = "MyData" # base of file name
Proteins = c()

#reading in complete RNA-seq data sets - has to be done only once!
RNASeqList = list()
SigDiffList = list()
for(PROTEIN in Proteins){
  InData=read.table(paste(FN,PROTEIN,"_vs_SCR.csv", sep = ""), sep = ",",header =T)
  RNASeqList[[PROTEIN]] = InData[order(InData$Row.names),]
  print(paste(PROTEIN,length(InData[,1])))
  SigDiffList[[PROTEIN]] = as.character(InData$Row.names[which(InData$padj <= 0.01)])
}
RNASeqDF = as.data.frame(do.call(cbind,RNASeqList))
RNASeqDF=RNASeqDF[c("X.Row.names", "X.external_gene_id","X.log2FoldChange","Y.log2FoldChange","Z.log2FoldChange")]
names(RNASeqDF)= c("GeneID","external_gene_id","X","Y","Z")

###################################################################################
#### run this once as well to get the function into your workspace

PlotMyAwesomeHeatmap <- function(KDs_of_interest= c("MSL1","MSL2","MOF"), Genes_of_interest = "test_ids_file.txt", I_want_a_pdf = TRUE, Out_File_Name = "MyAwesomeHeatmap",Only_Significantly_Affected = TRUE ){
  
  library(gplots)
  library(RColorBrewer)
  ids = read.table(Genes_of_interest)
  
  df = RNASeqDF[which(RNASeqDF$GeneID %in% ids$V1),] # if you want to extract based on sth else than the ENSMUS GeneID, change $GeneID to the correct column name, e.g. $external_gene_id
  
  if(Only_Significantly_Affected){
    for(PROTEIN in KDs_of_interest){
      nsa = which(!df$GeneID %in% SigDiffList[[PROTEIN]])
      for(i in nsa){
        df[i,][PROTEIN] = 0
      }
    }
  }
  
  ## this is to check for NAs and print out those that contain NA values
  containsNAs = df[!(complete.cases(df[KDs_of_interest])),]
  if(length(containsNAs[,1]) > 0){
    write.table(containsNAs,paste(Out_File_Name,"contained_NAs.txt", sep="_"), row.names = F, sep = "\t", quote = F)
    print("Check the file 'contained_NAs.txt' for genes that were omitted from the heatmap")
  }else{
    print("Congrats, no gene of interest contained NA's!")
  }
 # print("So far, you haven't generated a heatmap yet. Keep going and run the remaining lines.")
  
  ## generates matrix (without NA's) for heatmap
  forMM = df[complete.cases(df[KDs_of_interest]),]
  mm = as.matrix(forMM[KDs_of_interest])
  row.names(mm) = forMM$external_gene_id
  
  ## setting heatmap colors
  oranges <- colorRampPalette(brewer.pal(9,"Oranges"))(100)
  oranges <- rev(oranges)
  blues <- colorRampPalette(brewer.pal(9,"Blues"))(100)
  cols <- c(oranges,"grey", blues) ### if you don't like white for the color indicating the non-significantly changed genes, you can replace "white" by "black" or "grey" (or "green" or...)
  
  ## generates the heatmap
  if(I_want_a_pdf){
    pdf(paste(Out_File_Name,".pdf", sep =""))
    heatmap.2(mm, trace="none",margins=c(14,20), col=cols, symbreaks=TRUE)
    dev.off()
    print("Yay, you should have gotten a heatmap! Check your working directory!")
    } else{
      heatmap.2(mm, trace="none",margins=c(14,20), col=cols, symbreaks=TRUE)
    }
}

###################################################################################

# to generate a heatmap, type:
PlotMyAwesomeHeatmap(KDs_of_interest= c("MSL1","MSL2","MOF"), Genes_of_interest = "test_ids_file.txt", I_want_a_pdf = TRUE, Out_File_Name = "MyAwesomeHeatmap", Only_Significantly_Affected = TRUE)

# if you do not want to limit yourself to the values of significantly affected genes, the last option should say FALSE instead of TRUE
# if you want to see the plot immediatly in RStudio instead of getting a pdf, set I_want_a_pdf to FALSE
