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

colorSpectrum <- colorRampPalette(c("white",color))(100)
minValue=0
maxValue=max(metaData)/25
color_cuts <- seq(minValue,maxValue,length=101)

pdf(output)
heatmap.2(metaData, dendrogram='none', Rowv=FALSE, Colv=FALSE,
          trace='none', col=colorSpectrum, scale='none',
          labRow=c(), labCol=c(), breaks=color_cuts)
dev.off()


