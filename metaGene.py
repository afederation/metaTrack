import sys
sys.path.append('/ark/home/af661/src/utils/')
import utils

import os
import subprocess
import numpy


## For each locus, find the flanking region and divide it into 250 regions


def createRegionLoci(regionsLoci, extension, projectFolder, projectName):
    
    outputLoci = []
    locationToLocusDict = {}

    for locus in regionsLoci:
        
        location = locus.chr() + ':' + str(locus.start()) + '-' + str(locus.end())

        newlocus = utils.makeSearchLocus(locus, extension, extension)
        
        outputLoci.append(newlocus)
        locationToLocusDict[location] = locus

    return outputLoci, locationToLocusDict

def rankRegionsByAnchor(locationToLocusDict, anchorBam, numberRegions, projectFolder, projectName):
    
    print 'Ranking regions by anchor BAM'

    scoreDict = {}
    
    # Find the total number of reads in each region

    for location in locationToLocusDict:

        locus = locationToLocusDict[location]

        bamliquidatorCmd = 'bamliquidator %s %s %s %s . 1 200' % (anchorBam, locus.chr(),
                                                                  str(locus.start()), str(locus.end()))

        bamliquidatorOut = subprocess.Popen(bamliquidatorCmd, stdout = subprocess.PIPE, shell=True)
        score = bamliquidatorOut.communicate()[0]
        scoreDict[location] = int(score)

    # Rank the regions by the score and store in a dictionary
    
    allLocations = locationToLocusDict.keys()
    sortedLocations = sorted(allLocations, key=lambda x: scoreDict[x], reverse=True)
    
    rankDict = dict(zip(range(numberRegions), sortedLocations))
    
    return rankDict

def mapBamsToRegions(bamList, namesList, nbins, numberRegions, rankDict, locationToLocusDict, projectFolder, projectName):

    print 'Mapping BAMs to Regions'

    for b in range(len(bamList)):

        bam = bamList[b]
        sampleName = namesList[b]

        header = ['Region'] + range(nbins+1)[1:]
        data = [header]

        counter = 0

        for rank in range(numberRegions):
            
            counter += 1
            if counter%1000 == 0:
                print counter

            location = rankDict[rank]
            locus = locationToLocusDict[location]

            coords = (int(locus.start()), int(locus.end()))
            start = min(coords)
            stop = max(coords)

        
            cmd = ('bamliquidator ' + bam + ' ' + locus.chr() + ' ' + str(locus.start()) + ' '
                  + str(locus.end()) + ' . ' + str(nbins) + ' 0')

            bamliquidatorOut = subprocess.Popen(cmd, stdout = subprocess.PIPE, shell=True)
            rawOut = bamliquidatorOut.communicate()[0]
            binnedData = rawOut.split('\n')
            regionData = []

            for y in binnedData:
                if y != '':
                    regionData.append(int(y))
                else:
                    regionData.append(0)

            data.append(regionData)
    
        dataFilename = projectFolder + projectName + '_' + sampleName + '_metaData.txt'
        utils.unParseTable(data, dataFilename, '\t')

def makeGraphs(namesList, projectFolder, projectName):

    for sampleName in namesList:

        dataFilename = projectFolder + projectName + '_' + sampleName + '_metaData.txt'

        pdfFilename = projectFolder + projectName + '_' + sampleName + '_metaPlot.pdf'
        os.system('Rscript metaPlot.R ' + dataFilename + ' ' + pdfFilename)

        pdfFilename = projectFolder + projectName + '_' + sampleName + '_heatmap.pdf'
        os.system('Rscript heatmapPlot.R ' + dataFilename + ' ' + pdfFilename + ' red')



def main():
    
    from optparse import OptionParser

    usage = '''usage: %prog [options] -r [REGIONS_FILE] -b [BAM_FILES] -s [SAMPLE NAMES] 
             -a [ANCHOR DATASET] -o [OUTPUT DIRECTORY] -n [ANALYSIS NAME]'''
    parser = OptionParser(usage = usage)

    #required flags

    parser.add_option('-r','--regions', dest='regions',nargs = 1, default=None,
                      help = 'A file containing the regions of interest (BED or GFF)')
    parser.add_option('-b','--bams', dest='bams',nargs = 1, default=None,
                      help = 'Comma separated BAM files for analysis')
    parser.add_option('-s','--names', dest='names',nargs = 1, default=None,
                      help = 'Comma separated names that correspond to BAM files')
    parser.add_option('-a','--anchor', dest='anchor',nargs = 1, default=None,
                      help = 'The name used to rank regions by seq density')
    parser.add_option('-o','--output', dest='projectFolder',nargs = 1, default=None,
                      help = 'Output directory')
    parser.add_option('-n','--analysis_name', dest='projectName',nargs = 1, default=None,
                      help = 'Analysis Name')


    # Optional flags

    parser.add_option('','--bins', dest='bins',nargs = 1, default=250,
                      help = 'Number of bins for graphing regions')
    parser.add_option('','--extend', dest='extend',nargs = 1, default=1000,
                      help = 'Number of bases to extend the regions on both sides')
    parser.add_option('','--colors', dest='colors',nargs = 1, default='red',
                      help = 'List of colors to use, default is black')

    (options,args) = parser.parse_args()
    print(options)

    if (options.regions and options.bams and options.names and options.anchor
        and options.projectFolder and options.projectName):

        # Collect all the arguments

        regionsFile = options.regions
        regionsTable = utils.parseTable(regionsFile, '\t')
        regionsLoci = [utils.Locus(x[0], x[1], x[2], '.') for x in regionsTable]
        numberRegions = len(regionsLoci)

        projectFolder = options.projectFolder
        if projectFolder[-1] != '/':
            projectFolder += '/'
        projectName = options.projectName

        bams = options.bams
        bamList = bams.split(',')

        names = options.names
        namesList = names.split(',')

        anchor = options.anchor

        nbins = options.bins
        extension = options.extend
        colors = options.colors.split(',')

        try:
            anchorIndex = namesList.index(anchor)
        except ValueError:
            print 'ERROR: invalid anchor name'
        
        # Re-order the names and bams to put the anchor sample first

        anchorBam = bamList[anchorIndex]
        anchorName = namesList[anchorIndex]

        bamList.pop(anchorIndex)
        namesList.pop(anchorIndex)

        bamList.insert(0, anchorBam)
        namesList.insert(0, anchorName)

        # Run the functions

        extendedLoci, locationToLocusDict = createRegionLoci(regionsLoci, extension, projectFolder, projectName)
        
        rankDict = rankRegionsByAnchor(locationToLocusDict, anchorBam, numberRegions, projectFolder, projectName)

        mapBamsToRegions(bamList, namesList, nbins, numberRegions, rankDict, locationToLocusDict, projectFolder, projectName)

        makeGraphs(namesList, projectFolder, projectName)

    else:
        parser.print_help()
        sys.exit()



if __name__ == '__main__':
    main()
