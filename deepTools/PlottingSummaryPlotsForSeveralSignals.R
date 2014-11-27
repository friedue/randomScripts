## This script assumes that you have run computeMatrix at least twice with the --outFileNameData option using
## the _same_ BED-file, but different bigWig files.
## Please fill in the correct information in the next few lines 

#------------------------------------------------------
## required input: information about the data sets that were used with computeMatrix including labels and colors
#------------------------------------------------------

RegGroups = c("allGenes") # indicate the different bed-files that were used with computeMatrix, must be the same order as in the computeMatrix or profiler output, can be just one group, too
Signals= c("~/temp/dn57_H3K4me3_minArg_1000.tab","~/temp/wt26_H3K4me3_minArg_1000.tab") # list of paths to the files produced with computeMatrix (one for each bigWig, obtained with --outFileNameData), _must_ be more than one
SigLabels= c("dn","wt") # needs same order as the list of computeMatrix output files (Signals)
PlotTogether = "Signal" # either "Regions" or "Signal"; determines which two groups should be plotted within the same coordinate system (do you want to directly compare the coverage for different regions or different bigwig files); if you only had one group in the region file, it will not make sense to choose "Regions" here
GroupColors= c("green","red") # indicate the kind and the correct number of colors depending on how many different signals or regions you would like to plot within the same box

#------------------------------------------------------
## optional input
#------------------------------------------------------
PlotName = "MyAwesomeProfiles.png"
PlotTitle= "Summary Plots"
yAxisLabel = "bp"
xAxisLabel = "mean enrichment" 

## now run the lines above so that your specific info is read into R

#######################################################
## then run the following lines:
#######################################################
library(ggplot2) # you might need to install this library before you can use it (command: install.packages("ggplot2"))
source("~/temp/functions_PlottingSummaryPlotsForSeveralSignals.R")
PlottingSummaryPlots()
                                    
## you should now have a new .png file in your current working directory
