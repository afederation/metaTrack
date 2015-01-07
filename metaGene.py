import sys
sys.path.append('/ark/home/af661/src/utils/')
import utils

import os
import numpy

###
#
# Simplest case: finding the H3K27ac levels around p65 peaks
#
##

projectFolder = '/grail/projects/merkel/meta/'
projectName = 'MKL1'

## Bring in the regions of interest

regionsFile = '/grail/projects/merkel/macsEnriched/MKL1_p400Max_overlap.bed'
regionsTable = utils.parseTable(regionsFile, '\t')
regionsLoci = [utils.Locus(x[0], x[1], x[2], '.') for x in regionsTable]

## And the bams 

bamFile1 = '/grail/TONY/Decaprio/Homo.sapiens/Merkel.Cell.Carcinoma/MKL1/p400/20141229_2823/20141229_2823_hg19.sorted.bam'
bamFile2 = '/grail/TONY/Decaprio/Homo.sapiens/Merkel.Cell.Carcinoma/MKL1/Max/20141229_2824/20141229_2824_hg19.sorted.bam'

bamFiles = [bamFile1, bamFile2]
sampleNames = ['MAX', 'p400']
anchor = 'MAX'

## For each locus, find the flanking region and divide it into 250 regions

nbins = 250
extension = 1000

for b in range(len(bamFiles)):

    bam = bamFiles[b]
    sampleName = sampleNames[b]

    header = ['Region'] + range(nbins+1)[1:]
    data = [header]

    for locus in regionsLoci[1:100]:
    
        coords = (int(locus.start()), int(locus.end()))
        start = min(coords)
        stop = max(coords)

        print locus

        extendedLocus = utils.makeSearchLocus(locus, extension, extension)
        
        os.system('bamliquidator ' + bam + ' ' + extendedLocus.chr() + ' ' + str(extendedLocus.start()) + ' '
                  + str(extendedLocus.end()) + ' . ' + str(nbins) + ' 0 > ' + projectFolder + 'tempBamliquidator_'
                  + projectName + '.txt')    

        x = utils.parseTable(projectFolder + 'tempBamliquidator_' + projectName + '.txt', '\t')
        regionData = [int(y[0]) for y in x]

        regionLocation = [extendedLocus.chr() + ':' + str(extendedLocus.start()) + '-' + str(extendedLocus.end())]
        regionData = regionLocation + regionData

        data.append(regionData)

    sortedData = sorted(data[1:], key=lambda x: max([int(y) for y in x[1:]]), reverse=True)

    dataFilename = projectFolder + projectName + '_' + sampleName + '_metaData.txt'
    utils.unParseTable(sortedData, dataFilename, '\t')

    pdfFilename = projectFolder + projectName + '_' + sampleName + '_metaPlot.pdf'
    os.system('Rscript metaPlot.R ' + dataFilename + ' ' + pdfFilename)

    pdfFilename = projectFolder + projectName + '_' + sampleName + '_heatmap.pdf'
    os.system('Rscript heatmapPlot.R ' + dataFilename + ' ' + pdfFilename)
