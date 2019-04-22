#```{r exporting_tables, eval = FALSE}

#' Java garbage collection function
#' 
#' @description To avoid crashing of write.xlsx()
jgc <- function()
{
  gc()
  .jcall("java/lang/System", method = "gc")
} 

#--------------------------
## How to
#--------------------------

# java garbage collection
jgc()

### genes with peak numbers
library(xlsx)
write.xlsx(x = as.data.frame(genes.diffpeaks.out[comparison == "pg1_vs_cont",
                                                 -c("comparison"),
                                                 with =FALSE]),
           file = paste0("genes_with_differential_peaks_",
                         format(Sys.time(), "%Y%b%d"), ".xlsx"),
           sheetName = "pG1_vs_control", row.names = FALSE)

jgc()
