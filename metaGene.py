import sys
sys.path.append('/ark/home/af661/src/utils/')
import utils

import os
import subprocess
import numpy




def convertBEDtoGFF(regionsFile, extension, nbins, projectFolder, projectName):
    '''
    Take BED file from import
    Extend and divide each region given the paramaters
    Output a GFF with divided regions
    '''

    regionsTable = utils.parseTable(regionsFile, '\t')

    gff = []
    gffFilename = projectFolder + projectName + '_meta_regions.gff'
    
    counter = 0
    for line in regionsTable[:10000]:
        
        print counter
        counter += 1

        chrom = line[0]
        coords = (int(line[1]),  int(line[2]))
        start = min(coords)
        end = max(coords)

        # Check to see if additional information is present 
        try:
            name = line[3]
        except IndexError:
            name = projectName + '_' + str(counter)

        try:
            strand = line[5]
        except IndexError:
            strand = '.'


        regionLength = (end+extension) - (start-extension) 
        binLength = int(float(regionLength)/nbins)
        
        newStart = (start-extension)
        alternate = 1
        for n in range(nbins):

            ID = name + '|' + str(n)

            if alternate == 1:
                newEnd = newStart + binLength
            if alternate == -1:
                newEnd = newStart + binLength + 1


            gffLine = [chrom, ID, 'region', newStart, newEnd, 0, strand, '.', 1]

            newStart = newEnd
            alternate *= (-1)
            
            gff.append(gffLine)

    utils.unParseTable(gff, gffFilename, '\t')
    return gffFilename


def callBamliquidator(regionsGFFfile, projectFolder, namesList, bamList):

    for i in range(len(namesList)):

        name = namesList[i]
        bamFile = bamList[i]

        bamliquidatorFolder = projectFolder + name + '_liquidate/'
        utils.formatFolder(bamliquidatorFolder, True)

        mappingCmd = 'bamliquidator_batch'
        mappingCmd += ' -r ' + regionsGFFfile
        mappingCmd += ' -o ' + bamliquidatorFolder
        mappingCmd += ' -m -e 200 '
        mappingCmd += bamFile

        subprocess.call(mappingCmd, shell=True)
    
def parseBamliquidator(projectFolder, projectName, namesList, anchor, nbins):

    for name in namesList:

        bamliquidatorFile = projectFolder + name + '_liquidate/matrix.txt'
    
        parsedReads = []
        signalDict = {}
        dataDict = {}

        data = utils.parseTable(bamliquidatorFile, '\t')
           
        for line in data[1:]:
                
            IDstring = line[0].split('|')
            
            ID = IDstring[0]

            if ID not in signalDict:
                signalDict[ID] = 0
                dataDict[ID] = []
            signalDict[ID] += float(line[2])
            dataDict[ID].append(float(line[2]))

        if name == anchor:
            sortedID = sorted(signalDict.keys(), key=lambda x: signalDict[x], reverse=True)

        header = ['Region'] + range(nbins+1)[1:]
        output = [header]
        outputName = projectFolder + projectName + '_' + name + '_metaData.txt'

        for ID in sortedID:
            output.append([ID] + dataDict[ID])

        utils.unParseTable(output, outputName, '\t')


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


        projectFolder = options.projectFolder
        if projectFolder[-1] != '/':
            projectFolder += '/'
        projectName = options.projectName

        bams = options.bams
        bamList = bams.split(',')

        names = options.names
        namesList = names.split(',')

        anchor = options.anchor

        nbins = int(options.bins)
        extension = int(options.extend)
        colors = options.colors.split(',')


        anchorIndex = namesList.index(anchor)

        
        # Re-order the names and bams to put the anchor sample first

        anchorBam = bamList[anchorIndex]
        anchorName = namesList[anchorIndex]

        bamList.pop(anchorIndex)
        namesList.pop(anchorIndex)

        bamList.insert(0, anchorBam)
        namesList.insert(0, anchorName)

        # Run the functions

        regionsGFFfile = convertBEDtoGFF(regionsFile, extension, nbins, projectFolder, projectName)
        callBamliquidator(regionsGFFfile, projectFolder, namesList, bamList)
        parseBamliquidator(projectFolder, projectName, namesList, anchor, nbins)
        makeGraphs(namesList, projectFolder, projectName)

    else:
        parser.print_help()
        sys.exit()



if __name__ == '__main__':
    main()
