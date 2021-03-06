######
# 
# 150106
# MetaTrack R script 
# Plots data from corresponding python script
# Takes a set of regions and bamfiles, ranks and plots heatmaps
#     and meta-tracks of those regions
#
######

args = commandArgs(TRUE)
file = args[1]
output = args[2]
color = args[3]

library(gplots)


data = read.table(file)
metaData = as.matrix(data[,2:ncol(data)])
class(metaData) <- "numeric"

colorSpectrum <- colorRampPalette(c("white",color))(100)
minValue=0
maxValue=quantile(metaData,na.rm=TRUE,prob=0.5,names=FALSE)

color_cuts <- seq(minValue,maxValue,length=101)

png(output,width=1200,height=1200)
heatmap.2(metaData, dendrogram='none', Rowv=FALSE, Colv=FALSE,
          trace='none', col=colorSpectrum, scale='none',
          labRow=NA, labCol=NA, breaks=color_cuts)
dev.off()


