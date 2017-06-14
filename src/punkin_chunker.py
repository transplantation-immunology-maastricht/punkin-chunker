import sys
import os
from os.path import split, isfile, isdir, join, splitext
import subprocess
from subprocess import check_output
from os import makedirs
from Bio import SeqIO
from Blast_Result import Blast_Result
from Blast_Result_Set import Blast_Result_Set
from Bio.Blast.Applications import NcbiblastnCommandline
import operator


def CreateBlastDatabase(HLAReferenceFilename):
    print ('Creating a blast database.')
    makeBlastDB_cline = ('makeblastdb' 
        + ' -in ' + HLAReferenceFilename 
#        + ' -parse_seqids -dbtype nucl')
        + ' -dbtype nucl')
    print ('MakeDB Commandline:\n' + makeBlastDB_cline)
    subprocess.call(makeBlastDB_cline, shell=True)

# This method is a directory-safe way to open up a write file.
def createOutputFile(outputfileName):
    tempDir, tempFilename = split(outputfileName)
    if not isdir(tempDir):
        makedirs(tempDir)
    resultsOutput = open(outputfileName, 'w')
    return resultsOutput


def sortDirectory(readDirectory, outputDirectory, sortReference):
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
                print ('loading Reads from:\n' + fullReadFileName)
                #splitext[0] will get the filename without the .fastq extension.
                outputSubDirectory = join(outputDirectory,splitext(currentInputReadFile)[0])
                print('Sorted Reads will go here:\n' + outputSubDirectory)
                
                #BlastMinionReadsAgainstAPDRef(fullReadFileName, outputSubDirectory, sortReference)
                sortMinIONReadsByGene(fullReadFileName, outputSubDirectory, sortReference)
            
        else:
            # Nothing to do with this file for now
            pass    
        
            
            
            
            #newReads = createCollectionFromReadFile(join(inputReads,currentInputReadFile))
            #allReads.concatenate(newReads)
    

# TODO: I should detecte forward/reverse matches, and revcom th e sqeuence.  
# I have code to revcom the sequence in search_barcode, aka nit-picker now
def sortMinIONReadsByGene(inputFilename, outputDirectory, HLAReferenceFilename):
    print ('Time to sort our MinION Reads.  Compare against the APD Allele Reference.')
    
    # The full blast output can be several gigabytes worth of text.  Probably not worth writing.
    printFullBlastOutput = False

    # TODO : Method to read parameters, this is kind of crappy.  I know better.
    # I shouldn't be reading parameters at all.  
    #HLAReferenceFilename = sys.argv[1]
    #inputFilename = sys.argv[2]
    #outputDirectory = sys.argv[3]
    
    shortFilename = split(inputFilename)[1]
    #print('Short Filename:' + shortFilename)

    print ('HLA APD Reference:' + HLAReferenceFilename)
    print ('MinION Read Input file:' + inputFilename)
    print ('Output directory:' + outputDirectory)
    
    global FileOutputFormat
    if (".fasta" == inputFilename[-6:] or ".fa" == inputFilename[-3:]):
        FileOutputFormat = "fasta"
    elif (".fastq"== inputFilename[-6:] or ".fq" == inputFilename[-3:]):
        FileOutputFormat = "fastq"

    sortResultsOutput = createOutputFile(join(outputDirectory ,
        (shortFilename + '.SortResults.txt'))) 
    
    if(printFullBlastOutput):
        fullBlastOutput = createOutputFile(join(outputDirectory ,
            (shortFilename + '.FullBlastOutput.txt'))) 
    
    shortBlastOutput = createOutputFile(join(outputDirectory ,
        (shortFilename + '.ShortBlastOutput.txt'))) 
    
    finalBlastSummaryOutput = createOutputFile(join(outputDirectory ,
        (shortFilename + '.BlastSummary.txt'))) 
     
 
    print ('Parsing input file.  It\'s format is ' + FileOutputFormat)
    parsedReads = SeqIO.parse(inputFilename, FileOutputFormat)
    minionReadRecords = enumerate(parsedReads)
    readCount = len(list(SeqIO.parse(inputFilename, FileOutputFormat)))
    print (str(readCount) + ' reads found in input.')

    finalBlastSummaryOutput.write('HLA APD Reference:' + HLAReferenceFilename + '\n')
    finalBlastSummaryOutput.write('MinION Read Input file:' + inputFilename + '\n')
    finalBlastSummaryOutput.write('Output directory:' + outputDirectory + '\n')
    finalBlastSummaryOutput.write('Read Count: ' + str(readCount) + '\n')

    CreateBlastDatabase(HLAReferenceFilename)
    
    blastResultSets = []
 
    # Ineed to replace this for loop with 
    
    print ('Length of the enumerated read list:' + str(readCount))
 
    # TODO: I wonder if i can thread this blast stuff.  They all use the same database so probably not?
    # Maybe reads can be sorted in batches to be faster, because this takes forever.
    # Each record represents an HLA element in the input fasta file.
    for index, record in minionReadRecords:

        currentReadID = str(record.id)
        #print ('Read ID:' + currentReadID)
        currentSequence = str(record.seq)
        #print ('Read Sequence:' + currentSequence)

        print ('Sorting Read (' + str(index) + '/' + str(readCount) + ') : ' + currentReadID)
        
        currentBlastResultSet = Blast_Result_Set()
        currentBlastResultSet.readID = currentReadID
        currentBlastResultSet.readObject = record
        
        #print ('This fresh BlastResultSet should have 0 length blast results:' + str(len(currentBlastResultSet.BlastResults)))
        if(printFullBlastOutput):
            fullBlastOutput.write('\n\nBlasting Read:' + currentReadID + '\n')
        shortBlastOutput.write('\n\nBlasting Read:' + currentReadID + '\n')

        #Blast the read against the database
        # I can pass the sequence directly into the blast from stdio.  Pipe the sequence into blast.  
        # Make the commandline look like this:
        # echo -e ">Name\nGGTTGAATG" | blastn -outfmt 0 -db /home/ben/MUMCScripts/BlastMinIONReads/inputData/SimpleReference.fasta -evalue 0.001

        blastn_cline = NcbiblastnCommandline(db=HLAReferenceFilename, evalue=0.001, outfmt=0)
        commandLineQuery = 'echo -e \'>' + currentReadID + '\n' + currentSequence + '\' | ' + str(blastn_cline)

        # check_output will execute the commandline, wait for it to finish, and capture the output.
        blastStdIOText = check_output(commandLineQuery, shell=True)
        blastStdIOTextSplit = str(blastStdIOText).split('\n')
        
        blastResultList = list(blastStdIOTextSplit)
        blastResultLineCount = len(blastResultList)
        
        #print ('Blast Result Found, Line Count:' + str(blastResultLineCount))
            
        # Store Blast hits
        #for index, line in enumerate(blastStdIOTextSplit):
    
        blastResultLoopIndexer = 0
        while(blastResultLoopIndexer < blastResultLineCount):
        
            line = blastResultList[blastResultLoopIndexer] 


            #print(line)

            if(printFullBlastOutput):
                fullBlastOutput.write(line + '\n')
            
            if('>' in line):
                # The > character signifies this is the fasta header for the blast hit.  
                # This line contains the Allele name, which we can parse out.
                alleleNameLine = line
                shortBlastOutput.write(line)
                
                #print('AlleleNameFound:' + alleleNameLine)
                #print('Creating a brand new Blast Result.')
                
                currentBlastResult = Blast_Result()
                
                #print('Fresh Blast Result, should have 0 score ' + str(currentBlastResult.BlastScore))
                currentBlastResult.parseNomenclatureLine(alleleNameLine)
                
                foundScore = False

                while(not foundScore):
                    
                    blastResultLoopIndexer += 1
                    
                    line = blastResultList[blastResultLoopIndexer] 
                    
                    if(printFullBlastOutput):
                        fullBlastOutput.write(line + '\n')
                    
                    if('Score =' in line):
                        scoreLine = line
                        shortBlastOutput.write(line + '\n')
                        
                        #print('ScoreFound:' + scoreLine)
                        
                        #Assign the blast score
                        currentBlastResult.parseScoreLine(scoreLine)

                        foundScore = True
               
                currentBlastResultSet.BlastResults.append(currentBlastResult)  
                #print('Stored a blast result in my blast result set. Current count = ' + str(len(currentBlastResultSet.BlastResults)))     
                
            blastResultLoopIndexer += 1
        
        # Add the current blast result to this read's result set
        blastResultSets.append(currentBlastResultSet)
        #print('Stored a blast result set in my list of blast result sets. Current count = ' + str(len(blastResultSets)))



    #Loop through blastResultSets to print out some stats on each.
    for blastResultSet in blastResultSets:
        blastResultSet.assignReadToGeneAndAlleleGroup()
        blastResultSet.printResultSummary(sortResultsOutput)
        
        
    printSortedFastaFiles(blastResultSets, outputDirectory, finalBlastSummaryOutput)

    sortResultsOutput.close()
    if(printFullBlastOutput):
        fullBlastOutput.close()
    shortBlastOutput.close()
    
    finalBlastSummaryOutput.close()
    
    
# TODO: Maybe make some simple graphs to illustrate how reads were sorted.  Pie chart for Genes? Bar chart ?
def printSortedFastaFiles(ResultSets, outputDirectory, finalBlastSummaryOutput):
    sortedOutputDirectory = join(outputDirectory, 'SortedReads')
    # Store Tuples: (Gene, OutputFile, readCount)
    geneLevelOutputFiles = []
    # Store Tuples: (Gene, Group, OutputFile, readCount)
    groupLevelOutputFiles = []
    
    unsortedReadOutputFile = createOutputFile(join(sortedOutputDirectory,'UnsortedReads.' + FileOutputFormat))
    unsortedReadCount = 0

    global FileOutputFormat
    
    for currentResultSet in ResultSets:
        currentGene = currentResultSet.AssignedGene
        currentGroup = currentResultSet.AssignedAlleleGroup
        
        
        #Here's the logic I'm about to accomplish:
        
        # If we have the gene assigned:
            # If we have the group assigned:
                # Write file to group level output.
            # Else
                # Write file to Gene level output.
        # Else
            # Write to rejected read output.
        
        
        
           
        # If the gene is assigned   
        if(len(currentGene) > 0 and currentGene != '-1'):
        
        

            # If the allele Group (first field) is assigned
            if(len(currentGroup) > 0 and currentGroup != '-1'):
        
                currentGroupLevelOutputFile = None
                
                # Search for existing Gene level outputffiles.
                foundGroupLevelOutput = False
                
                # If we already have a group level output file in our list, use that one.
                for groupLevelIndex,groupLevelOutputFile in enumerate(groupLevelOutputFiles):
                    if (groupLevelOutputFile[0] == currentGene and groupLevelOutputFile[1] == currentGroup):
                        foundGroupLevelOutput = True
                        currentGroupLevelOutputFile = groupLevelOutputFile[2]
                        
                        # This will increment the read count, by replacing the whole groupLevelOutputFile tuple.
                        # There is certainly a better way to do this.
                        groupLevelOutputFiles[groupLevelIndex] = (
                            groupLevelOutputFile[0],
                            groupLevelOutputFile[1],
                            groupLevelOutputFile[2],                            
                            groupLevelOutputFile[3] + 1)
                
                # None found, make a new output file.
                if not foundGroupLevelOutput:
                    currentGroupLevelOutputFile = createOutputFile(join(sortedOutputDirectory,'HLA-' + currentGene + '_' + currentGroup + '.' + FileOutputFormat))
                    groupLevelOutputFiles.append((currentGene,currentGroup,currentGroupLevelOutputFile,1))
                        
                # Print the sequence to the Group level output.
                if foundGroupLevelOutput != None:
                    SeqIO.write([currentResultSet.readObject], currentGroupLevelOutputFile, FileOutputFormat)
                    
            else: 
                print ('This read maps to a gene but not a group:' + str(currentResultSet.readID))
                
                # Print the read to a group level             
                currentGeneLevelOutputFile = None
                
                # Search for existing Gene level outputffiles.
                foundGeneLevelOutput = False
                
                # If we already have a gene level output file in our list, use that one.
                for geneLevelIndex, geneLevelOutputFile in enumerate(geneLevelOutputFiles):
                    if geneLevelOutputFile[0] == currentGene:
                        foundGeneLevelOutput = True
                        currentGeneLevelOutputFile = geneLevelOutputFile[1]
                        
                        geneLevelOutputFiles[geneLevelIndex] = (
                                geneLevelOutputFile[0],
                                geneLevelOutputFile[1],
                                geneLevelOutputFile[2] + 1)
                        
                # None found, make a new output file.
                if not foundGeneLevelOutput:
                    currentGeneLevelOutputFile = createOutputFile(join(sortedOutputDirectory,'HLA-' + currentGene + '.' + FileOutputFormat))
                    geneLevelOutputFiles.append((currentGene,currentGeneLevelOutputFile,1))
                        
                # Print the sequence to the Gene level output.
                if foundGeneLevelOutput != None:
                    SeqIO.write([currentResultSet.readObject], currentGeneLevelOutputFile, FileOutputFormat)
                
                    
        # This else corresponds to if there was no gene assigned for this read.  
        # I need to do something with the unsorted reads.  
        else:
            unsortedReadCount += 1
            SeqIO.write([currentResultSet.readObject], unsortedReadOutputFile, FileOutputFormat)            
            #print ('YOU NEED TO DO SOMETHING WITH THIS READ WHICH DOES NOT MAP TO THE REFERENCE:' + str(currentResultSet.readID))
            
        
    # Write sort results to the output file and close the Gene and Group specific fasta files.    
    finalBlastSummaryOutput.write('\n\nSorting Read Results:\n')

    #for geneLevelOutputFile in geneLevelOutputFiles:
    for geneLevelOutputFile in sorted(geneLevelOutputFiles, key=operator.itemgetter(2), reverse=True): 
        finalBlastSummaryOutput.write('HLA-' + geneLevelOutputFile[0] + ': ' + str(geneLevelOutputFile[2]) + ' Reads\n')
        geneLevelOutputFile[1].close()
        
    #for groupLevelOutputFile in groupLevelOutputFiles:
    #Sorted by number of reads
    for groupLevelOutputFile in sorted(groupLevelOutputFiles, key=operator.itemgetter(3), reverse=True): 
        finalBlastSummaryOutput.write('HLA-' + groupLevelOutputFile[0] +  groupLevelOutputFile[1] + ': ' + str(groupLevelOutputFile[3]) + ' Reads\n')
        groupLevelOutputFile[2].close()
    
    
    finalBlastSummaryOutput.write('Unsorted Reads: ' + str(unsortedReadCount) + ' Reads\n')
    unsortedReadOutputFile.close()
    

