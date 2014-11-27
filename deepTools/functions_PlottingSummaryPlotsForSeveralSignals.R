ReadingInFidelsSummaryPlots = function(File, Groupsname, Groups, Signalname){
          InData = read.table(File)
          GroupData_List = list()
          k = 2
          for(i in c(1:length(Groups))){
            GroupData = InData[c(1,k,k+1)]
            names(GroupData)= c("bin","mean","std")
            GroupData[Groupsname] = Groups[i]
            k = k + 2
            GroupData_List[[i]] = GroupData
          }
          
          AllGroups = do.call(rbind, GroupData_List)
          AllGroups$Signal = Signalname
          
          return(AllGroups)
        }
        
PlottingSummaryPlots = function(ListOfInputFiles =  Signals, RegionGroups = RegGroups, Signalnames = SigLabels, Group = PlotTogether,
                                    OutPlotName = PlotName, Title = PlotTitle, YLab = yAxisLabel, XLab = xAxisLabel ){
    SumPlots=list()
    i=1
    for(s in ListOfInputFiles){
       InData = ReadingInFidelsSummaryPlots(s,"Regions",RegionGroups,Signalnames[i])
       SumPlots[[i]] = InData
       i=i+1
        }
    
    SumPlots_df = as.data.frame(do.call(rbind, SumPlots))
    
if(Group=="Regions"){
    PLOT = ggplot() +  geom_line(data = SumPlots_df, aes(x=bin, y = mean, colour = Regions)) +
            theme_bw(base_size = 16) +
            facet_grid(Signal~.) +
            scale_colour_manual(values=GroupColors,name="") +
            theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
             ylab(YLab) + xlab(XLab)+
            ggtitle(Title)
            }
else{PLOT = ggplot() +  geom_line(data = SumPlots_df, aes(x=bin, y = mean, colour = Signal)) +
            theme_bw(base_size = 16) +
            facet_grid(Regions~.) +
            scale_colour_manual(values=GroupColors,name="") +
            theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
            ylab(YLab) + xlab(XLab)+
            ggtitle(Title)
            }
  
    png(OutPlotName) # change this if you want a pdf, tiff etc.
    print(PLOT)
    devname=dev.off()
}
