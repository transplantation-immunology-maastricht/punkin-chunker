import sys
import os
from os.path import split, isfile, isdir, join, splitext
import subprocess
#from subprocess import check_output
from subprocess import Popen, PIPE, STDOUT
from os import makedirs
from Bio import SeqIO
from Blast_Result import Blast_Result
#from Blast_Result_Set import Blast_Result_Set
#from Bio.Blast.Applications import NcbiblastnCommandline
import operator
from Bio.Blast import NCBIXML
import StringIO
import datetime


def CreateBlastDatabase(HLAReferenceFilename):
    #print ('Creating a blast database...')
    makeBlastDB_cline = ('makeblastdb' 
        + ' -in ' + HLAReferenceFilename 
#        + ' -parse_seqids -dbtype nucl')
        + ' -dbtype nucl')
    #print ('MakeDB Commandline:\n' + makeBlastDB_cline)
    subprocess.call(makeBlastDB_cline, shell=True)

# This method is a directory-safe way to open up a write file.
def createOutputFile(outputfileName):
    tempDir, tempFilename = split(outputfileName)
    if not isdir(tempDir):
        makedirs(tempDir)
    resultsOutput = open(outputfileName, 'w')
    return resultsOutput


def sortDirectory(readDirectory, outputDirectory, sortReference, threadCount):
    print ('Sorting this directory of reads:\n' + str(readDirectory))
    
    for currentInputReadFile in os.listdir(readDirectory):
        readFormat = None
        if (".fasta" == currentInputReadFile[-6:] or ".fa" == currentInputReadFile[-3:]):
            readFormat = 'fasta'
            # TODO Fix this, support fasta.
            raise Exception('Bad Read Input Format, I dont support fasta in the blast sorter yet')
        elif(".fastq"== currentInputReadFile[-6:] or ".fq" == currentInputReadFile[-3:]):
            readFormat = 'fastq'
        else :
            print ('Skipping this file:' + currentInputReadFile)
            #raise Exception('Bad Read Input Format')
        
        if(readFormat is not None):
            # Skip files that contain '_Rejected_'
            if('_Rejected_' in currentInputReadFile):
                print ('Skipping Rejected Read File:\n' + currentInputReadFile)
            elif('_Unfiltered_' in currentInputReadFile):
                print('Skipping Unfiltered Read File:\n' + currentInputReadFile)
            else:
                fullReadFileName=str(join(readDirectory,currentInputReadFile))
                #print ('loading Reads from:\n' + fullReadFileName)
                #splitext[0] will get the filename without the .fastq extension.
                outputSubDirectory = join(outputDirectory,splitext(currentInputReadFile)[0])
                print('Sorted Reads will go here:\n' + outputSubDirectory)
                
                #BlastMinionReadsAgainstAPDRef(fullReadFileName, outputSubDirectory, sortReference)
                sortMinIONReadsByGene(fullReadFileName, outputSubDirectory, sortReference, threadCount)
            
        else:
            # Nothing to do with this file for now
            pass    

 


def sortMinIONReadsByGene(inputFilename, outputDirectory, HLAReferenceFilename, threadCount):
    #print ('Time to sort our MinION Reads.  Compare against the APD Allele Reference.')

    shortFilename = split(inputFilename)[1]

    print ('HLA APD Reference:\n' + HLAReferenceFilename)
    print ('MinION Read Input file:\n' + inputFilename)
    print ('Output directory:\n' + outputDirectory)
    
    global FileOutputFormat
    if (".fasta" == inputFilename[-6:] or ".fa" == inputFilename[-3:]):
        FileOutputFormat = "fasta"
    elif (".fastq"== inputFilename[-6:] or ".fq" == inputFilename[-3:]):
        FileOutputFormat = "fastq"

    #sortResultsOutput = createOutputFile(join(outputDirectory ,
    #    (shortFilename + '.SortResults.txt'))) 
    
    #shortBlastOutput = createOutputFile(join(outputDirectory ,
    #    (shortFilename + '.ShortBlastOutput.txt'))) 
    
    finalBlastSummaryOutput = createOutputFile(join(outputDirectory ,
        (shortFilename + '.BlastSummary.txt'))) 
     
 
    #print ('Parsing input file.  It\'s format is ' + FileOutputFormat)
    parsedReads = SeqIO.parse(inputFilename, FileOutputFormat)
    #minionReadRecords = enumerate(parsedReads)
    #readRecordList = list(minionReadRecords)
    readRecordList = list(parsedReads)

    readCount = len(readRecordList)
    print (str(readCount) + ' reads ready to sort.')
    
    #for i in range(0,10):
    #    print ('index ' + str(i) + 'readrecord at this position is:' + str(readRecordList[i]))
  
    # pause statement for debugging...
    #wait = raw_input("\n\n\n\n**************************PRESS ENTER TO CONTINUE.**************************\n\n\n\n")

    # Print some info to the sorting summary
    finalBlastSummaryOutput.write('HLA APD Reference:' + HLAReferenceFilename + '\n')
    finalBlastSummaryOutput.write('MinION Read Input file:' + inputFilename + '\n')
    finalBlastSummaryOutput.write('Output directory:' + outputDirectory + '\n')
    finalBlastSummaryOutput.write('Read Count: ' + str(readCount) + '\n')

    # Create the Blast Reference before we start.
    CreateBlastDatabase(HLAReferenceFilename)


    # Split the reads into batches
    batchSize = 500
    
    readBatches = splitReadsIntoBatches(readRecordList, batchSize)
    # After splitting I'll clear the old value, because memory management. This is probably unnecessary, compilers are smart.
    readRecordList = None
    
    print ('I split the ' + str(readCount) + ' reads into ' 
        + str(len(readBatches)) + ' batches of <=' + str(batchSize) + ' reads.')
    

    # pause statement for debugging...
    #wait = raw_input("\n\n\n\n**************************PRESS ENTER TO CONTINUE.**************************\n\n\n\n")

    
    combinedResultSets = []
    
    # This loop processes each batch of reads indivudually
    # TODO: Add some threading to this.
    for batchIndex, readBatch in enumerate(readBatches):
        print 'Sorting Batch # (' + str(batchIndex + 1) + '/' + str(len(readBatches)) + ') ' + str(datetime.datetime.now())
        
        currentResults = sortReadArrayByGene(readBatch, HLAReferenceFilename, threadCount)
        
        #print ('I received ' + str(len(currentResults)) + ' results from batch # (' + str(batchIndex + 1) + '/' + str(len(readBatches)) + ')') 
        combinedResultSets = combinedResultSets + currentResults

    #blastResultSets = sortReadArrayByGene(readRecordList, HLAReferenceFilename, sortResultsOutput, shortBlastOutput)
    
    writeSortedReadsToFile(combinedResultSets, outputDirectory, finalBlastSummaryOutput)
    
        
        
    #printSortedFastaFiles(blastResultSets, outputDirectory, finalBlastSummaryOutput)

    #sortResultsOutput.close()
    #shortBlastOutput.close()
    
    finalBlastSummaryOutput.close()
    


def splitReadsIntoBatches(recordList, batchSize):
    #print ('I have ' + str(len(recordList)) + ' reads to split into batches of size ' + str(batchSize))
    # remainderRecords keeps track of the "unbatched" reads.
    remainderRecords = recordList
    # readBatches is an array full of read-arrays.  Sort of a 2d array but not quite.
    readBatches = []
        
    while len(remainderRecords) > 0:
        currentBatch = remainderRecords[0:batchSize]
        readBatches.append(currentBatch)
        remainderRecords = remainderRecords[batchSize:]

    return readBatches

def sortReadArrayByGene(minionReadRecordList, HLAReferenceFilename, threadCount):
    # TODO: I wonder if i can thread this blast stuff.  They all use the same database so probably not?
    # Maybe reads can be sorted in batches to be faster, because this takes forever.
    # Each minionReadRecord represents an HLA element in the input fasta file.
    # TODO: Now i think, instead of threading, I should Batch the reads into groups, and tell blastn to thread itself.
    # I bet blastn is better at threading than I am.
    
    
    # I will move this a layer higher to be more efficient.
    #CreateBlastDatabase(HLAReferenceFilename)
    
    #readCount = len(minionReadRecordList)
    
    blastResultSet = []

    # I group the reads into a single fasta string to pass to BLAST  
    fastaString = ''
        
    for minionReadIndex, minionReadRecord in enumerate(minionReadRecordList):
        currentReadID = str(minionReadRecord.id)
        currentSequence = str(minionReadRecord.seq) 
        
        fastaString += '>' + currentReadID + '\n' + currentSequence + '\n'

    # TODO: this thread count is not really working.  
    # TODO: I should do real threading by myself.  
    commandLineArray = ['blastn'
        , '-db' , HLAReferenceFilename               
        , '-num_threads', str(threadCount)         
        , '-outfmt', '5'
        , '-evalue', '0.001'             
        ]
    
    #print('COMMANDLINEARRAY')
    #print(str(commandLineArray))
    
    blastProcessHandle = Popen(commandLineArray, stdout=PIPE, stdin=PIPE, stderr=STDOUT)    
    #grep_stdout = p.communicate(input=fastaString)[0]

    blastXMLText = blastProcessHandle.communicate(input=fastaString)[0].decode()
    
    #print ('resulting blast xml text:')
    #print(blastXMLText)
    
    blastResults = parseXMLForBlastResults(minionReadRecordList,blastXMLText)
    
    # I don't use append, because, blastResults might have more than one minionReadRecord.  Might.
    blastResultSet = blastResultSet + (blastResults)
        
        
    #print('I am done, and I need to do a final check.')
    #print('Length MinION Record List: ' + str(len(minionReadRecordList)))
    #print('Length MblastResultSet: ' + str(len(blastResultSet)))



    return blastResultSet 

# A method to print information about blast alignments.
# I think this method is unneccsary. It just prints alignment stuff from
# the blast xml

           

def parseXMLForBlastResults(readRecords, blastXMLText):
    
    #print('parsing xml:')
    #print('Length readREcords =' + str(len(readRecords)))
    
    
    # pause statement for debugging...
    #wait = raw_input("\n\n\n\n**************************PRESS ENTER TO CONTINUE.**************************\n\n\n\n")

    #Parse XML for blast results. Need a File handle first.
    xmlStringHandle = StringIO.StringIO(blastXMLText)    
    blastRecords = NCBIXML.parse(xmlStringHandle)
    
    #print ('After Blast:')
    #print ('Length readRecords:' + str(len(readRecords)))
    #print ('Length blastRecords:' + str(len(blastRecords)))

    blastResults = []
    
    # parse returns multiple blast results.
    # Right now Im blasting one read at a time.  
    # I can use parse if I somehow feed multiple
    # reads to BLAST (which I should, I'm sure blast is efficienter than me))
    for blastRecordIndex, blastRecord in enumerate(blastRecords):  
        
        currentBlastResult = Blast_Result()   
        currentBlastResult.readRecord = readRecords[blastRecordIndex]
        
        if(len(blastRecord.alignments) < 1):
            #print ('No alignments detected for this read. This is not a problem.')
            pass
        else:
            currentBlastResult.processBlastRecord(blastRecord)
                    
        blastResults.append(currentBlastResult)  


    return blastResults



#I think this method is deprecated.
def printBlastRecordInformation(blastRecord):
    
    print ('BlastRecord:' + str(blastRecord))
    #print ('this record aligned with: ' + blast)
    
    currentAlignments = blastRecord.alignments
    print ('Number of alignments:' + str(len(currentAlignments)))
    for alignment in currentAlignments:
        print '****Alignment****'
        print 'sequence:' + str(alignment.title)
        print 'length:' + str(alignment.length)
        # Usually there is only one HSP.
        # hsp stores info about the alignment
        currentHsps = alignment.hsps
        if (len(currentHsps) != 1):
            print 'I only expect exactly one hsp for this alignment. There are ' + str()
            raise Exception('Multiple HSPs. I don\'t know what this means but you have to fix it.')
        #print ('Number of hsps:' + str(len(currentHsps)))
        # I don't need to loop these. Just grab the first.
        for hsp in currentHsps: #print('e value:' + str(hsp.expect))
            print 'score:' + str(hsp.score)
            print 'strand:' + str(hsp.strand)
            print 'frame:' + str(hsp.frame) #print('strand:' + str(hsp.strand))
            print hsp.query[0:75] + '...'
            print hsp.match[0:75] + '...'
            print hsp.sbjct[0:75] + '...'
            
# readRecordList contains a list of the reads.
# blastResultSets is a list of the blastResults.
# TODO: I really should sync those together.  Maybe put the read record right into the blast result object.  Yeah do that.
# Because It's a pain to keep the two arrays in sync.   With threading it's much harder.
# I Think i did this, but i'll have to watch to see what i broke.
def writeSortedReadsToFile(blastResults, outputDirectory, finalBlastSummaryOutput):    
    print('I will now write the sorted reads to file.')
    
    #sortedOutputDirectory = join(outputDirectory, 'SortedReads')
    
    # Store Tuples: (Gene, OutputFile, readCount)
    geneLevelOutputFiles = []

    unsortedReadOutputFile = createOutputFile(join(outputDirectory,'UnsortedReads.' + FileOutputFormat))
    unsortedReadCount = 0

    global FileOutputFormat

    for resultIndex, currentBlastResult in enumerate(blastResults):
        
        if(currentBlastResult is None):
            # This corresponds to if there was no blast results for this read. 
            unsortedReadCount += 1
            SeqIO.write([currentBlastResult.readRecord], unsortedReadOutputFile, FileOutputFormat)
            raise Exception('This read was unsorted.  No problem, but I wanted to see when this happens. Delete this exception.')
               
        else:
            currentGene = currentBlastResult.Gene
            #currentGroup = currentBlastResult.AssignedAlleleGroup
    
            # If the gene is assigned   
            if(currentGene is not None and len(currentGene) > 0 and currentGene != '-1'):
         
                # Print the read to a group level             
                currentGeneLevelOutputFile = None
                
                # Search for existing Gene level outputffiles.
                #foundGeneLevelOutput = False
                
                # If we already have a gene level output file in our list, use that one.
                # I keep track of the outputFileIndex on my own, because enumerate() is smarter than I am.
                outputFileIndex = 0
                for outputFileGene, outputFileObject, readCount in geneLevelOutputFiles:
                    #print ('outputFileObject is this:' + str(outputFileObject))
                    #print ('outputFileObject[0] is this:' + str(outputFileObject[0]))
                    if outputFileGene == currentGene:
                        #foundGeneLevelOutput = True
                        currentGeneLevelOutputFile = outputFileObject
                        
                        # Increment Read Count
                        geneLevelOutputFiles[outputFileIndex] = (
                            outputFileGene,
                            outputFileObject,
                            readCount + 1)
                    
                    outputFileIndex += 1
                # None found, make a new output file.
                if currentGeneLevelOutputFile is None:
                    currentGeneLevelOutputFile = createOutputFile(join(outputDirectory,'HLA-' + currentGene + '.' + FileOutputFormat))
                    geneLevelOutputFiles.append((currentGene,currentGeneLevelOutputFile,1))
                        
                # Print the sequence to the Gene level output.
                if currentGeneLevelOutputFile != None:
                    # Is the sequence mapped to the reference in the reverse direction?
                    
                    if (currentBlastResult.ForwardMatch):
                        SeqIO.write([currentBlastResult.readRecord], currentGeneLevelOutputFile, FileOutputFormat)
                    else:
                        #print ('about to check the reverse complement of this read:\n')
                        #print (str(currentBlastResult.readRecord))
                        reverseRecord = currentBlastResult.readRecord
                        forwardRecord = reverseRecord.reverse_complement(id=reverseRecord.id+"_reverse_complement", name=True, description=True)
                        #SeqIO.write([rec], allReadFileOutputs[barcodeKey], 'fastq')
                
                        #print ('This record is in the reverse direction.')
                        SeqIO.write([forwardRecord], currentGeneLevelOutputFile, FileOutputFormat)
                    
     
        
    # Write sort results to the output file and close the Gene and Group specific fasta files.    
    finalBlastSummaryOutput.write('\n\nSorting Read Results:\n')

    #for outputFileObject in geneLevelOutputFiles:
    for outputFileObject in sorted(geneLevelOutputFiles, key=operator.itemgetter(2), reverse=True): 
        finalBlastSummaryOutput.write('HLA-' + outputFileObject[0] + ': ' + str(outputFileObject[2]) + ' Reads\n')
        outputFileObject[1].close()

    finalBlastSummaryOutput.write('Unsorted Reads: ' + str(unsortedReadCount) + ' Reads\n')
    unsortedReadOutputFile.close()

    

