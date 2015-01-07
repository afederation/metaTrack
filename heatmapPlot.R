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

data = read.table(file)
metaData = as.matrix(data[,2:ncol(data)])

colorSpectrum <- colorRampPalette(c("white","red"))(100)
minValue <- quantile(metaData,na.rm=TRUE,prob=0.3,names=FALSE)
maxValue <- quantile(metaData,na.rm=TRUE,prob=0.90,names=FALSE)
color_cuts <- seq(minValue,maxValue,length=100)
color_cuts <- c(0, color_cuts,max(metaData))

pdf(output)
heatmap(metaData, Rowv=NA, Colv=NA, labRow=NA, labCol=NA, col=colorSpectrum)
dev.off()