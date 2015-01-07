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

metaData = colMeans(data[,2:ncol(data)])

pdf(output)

plot(metaData, type='l', xlab='bin', ylab='rpm/bp')

dev.off()