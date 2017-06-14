import operator

#acceptableScoreCutoff = 650
#acceptableScoreCutoff = 375
acceptableScoreCutoff = 100

class Blast_Result_Set:
    readID = ''
    readObject = None
    BlastResults = []
    AlleleGroupMatches = []
    AssignedGene = ''
    AssignedAlleleGroup = ''
    
    def __init__(self):
        self.readID = ''
        self.readObject = None
        self.BlastResults = []
        self.AlleleGroupMatches = []
        self.AssignedGene = ''
        self.AssignedAlleleGroup = ''
        
#    def assignReadToGene(self):
#        self.AssignedGene = 'A'


    
    def assignReadToGeneAndAlleleGroup(self):
        
        self.AssignedGene = '-1'
        self.AssignedAlleleGroup = '-1'
        foundGenes = self.uniqueGenes()       
        
        #print('Length of Found Genes:' + str(len(foundGenes)))
        
        if(len(foundGenes) != 1):
            # No match found.  
            #print('No Gene found for read:' + self.readID)

            pass
        
        
        
        else:
            # One Gene found, we can probably assign the gene confidently.
            #print('Assigning gene ' + str(foundGenes[0][0]) + ' for read:' + self.readID)
            self.AssignedGene = str(foundGenes[0][0])
            
            foundAlleleGroups = self.uniqueAlleleGroups(str(foundGenes[0][0]))  
            
            maxGroupScore = foundAlleleGroups[0][2]
            #print('Max Score found ' + str(maxGroupScore) + ' for read:' + self.readID)
            
            if(maxGroupScore > acceptableScoreCutoff):
                self.AssignedAlleleGroup = str(foundAlleleGroups[0][0])
            
            
            
        
        #for gene in foundGenes:
        #    foundAlleleGroups = self.uniqueAlleleGroups(str(gene[0]))            
        #    for alleleGroup in foundAlleleGroups:
        #        outputFile.write('\t\tHLA-' + str(gene[0]) + ':' + str(alleleGroup[0])  + ' : ' + str(alleleGroup[1])
        #            + ' hits,\tAverage Score=' + str(alleleGroup[2]) + '\n')   
                
        
        #self.AssignedGene = 'A'
        #self.AssignedAlleleGroup = '01'
    
    def printResultSummary(self, outputFile):
        outputFile.write('\n' + self.readID + '\n')
        
        if(self.AssignedGene == '' or self.AssignedGene == '-1'):
            outputFile.write('No Gene Assigned\n')
        else:
            outputFile.write('Assigned Gene ' + self.AssignedGene + '\n')
            
        if(self.AssignedAlleleGroup == '' or self.AssignedAlleleGroup == '-1'):
            outputFile.write('No Group Assigned\n')
        else:
            outputFile.write('Assigned Group ' + self.AssignedAlleleGroup + '\n')
        
        #print('Summarizing results for read ' + self.readID)
        
        foundGenes = self.uniqueGenes()
        
        
        for gene in foundGenes:
            #print('GeneFound:' + gene)
            #outputFile.write('\tHLA-' + gene + ' : ' + str(self.countBlastHitsByGene(gene)) 
            #    + ' hits,\tAverage Score=' + str(self.averageScoreByGene(gene)) + '\n')
            outputFile.write('\tHLA-' + str(gene[0]) + ' : ' + str(gene[1]) 
                + ' hits,\tAverage Score=' + str(gene[2]) + '\n')
            
            foundAlleleGroups = self.uniqueAlleleGroups(str(gene[0]))
            
            for alleleGroup in foundAlleleGroups:
                outputFile.write('\t\tHLA-' + str(gene[0]) + ':' + str(alleleGroup[0])  + ' : ' + str(alleleGroup[1])
                    + ' hits,\tAverage Score=' + str(alleleGroup[2]) + '\n')   
                #outputFile.write('\t\tHLA-' + gene + ':' + alleleGroup  + ' : ' + str(self.countBlastHitsByGroup(gene, alleleGroup)) 
                #    + ' hits,\tAverage Score=' + str(self.averageScoreByGroup(gene, alleleGroup)) + '\n')   
                #outputFile.write('\t\tHLA-' + gene + ':' + alleleGroup  + ' : ' + str(self.countBlastHitsByGroup(gene, alleleGroup)) 
                #    + ' hits,\tAverage Score=' + str(foundAlleleGroups[alleleGroup]) + '\n')     
                
                   
    def shortBlastResults(self):
        shortBlastCount = 9999999
        if(len(self.BlastResults) < shortBlastCount):
            return  self.BlastResults
        else:
            #print('returning shorter blast count.')
            return self.BlastResults[0:shortBlastCount]
    
    def uniqueGenes(self):
        # Store it in a list of tuples.  
        # Gene:NumberRecords:AvgBlastScore
        uniqueGeneNames = []
        
        #print('length of blast results:' + str(len(self.BlastResults)))
        #print('length of short blast results:' + str(len(self.shortBlastResults())))
        
        for BlastResult in self.shortBlastResults():
            #print('inside loop?')
            #if BlastResult.Gene not in uniqueGeneNames:
            if not filter(lambda a: a[0] == BlastResult.Gene, uniqueGeneNames):    
                #print('new gene found:' + BlastResult.Gene)
                #uniqueGeneNames.append(BlastResult.Gene)
                #uniqueGeneNames[BlastResult.Gene] = self.averageScoreByGene(BlastResult.Gene)
                uniqueGeneNames.append((BlastResult.Gene, self.countBlastHitsByGene(BlastResult.Gene) , self.averageScoreByGene(BlastResult.Gene)))
            
            #else:
                #print('this gene is not new:' + BlastResult.Gene)
            
        #print('Number Unique Genes:' + str(len(uniqueGeneNames)))
        # Sort by average Blast score descending.
        #print('Before sorting, length is' + str(len(uniqueGeneNames)))
        sortedList = sorted(uniqueGeneNames, key=operator.itemgetter(2), reverse=True)
        #print('After sorting, length is' + str(len(sortedList)))
        #return uniqueAlleleGroupNames
        return list(sortedList)
        
    def uniqueAlleleGroups(self,Gene):
        # Store it in a list of tuples.  
        # AlleleGroup:NumberRecords:AvgBlastScore
        uniqueAlleleGroupNames = []
        
        for BlastResult in self.shortBlastResults():
            if BlastResult.Gene == Gene:
                #uniqueAlleleGroupNames contains a tuple of data.  I need to fix the check for existence.
                if not filter(lambda a: a[0] == BlastResult.NomenclatureFields[0], uniqueAlleleGroupNames):             
                #if BlastResult.NomenclatureFields[0] not in uniqueAlleleGroupNames:
                    # Store Allele Group and Avg Score
                    #uniqueAlleleGroupNames[BlastResult.NomenclatureFields[0]] = self.averageScoreByGroup(Gene, BlastResult.NomenclatureFields[0])
                    uniqueAlleleGroupNames.append((BlastResult.NomenclatureFields[0] , self.countBlastHitsByGroup(Gene, BlastResult.NomenclatureFields[0]) , self.averageScoreByGroup(Gene, BlastResult.NomenclatureFields[0])))
            
        # Sort by average Blast score descending.
        sortedList = sorted(uniqueAlleleGroupNames, key=operator.itemgetter(2), reverse=True)
        return list(sortedList)
        
    def countBlastHitsByGene(self, gene):
        hitCounter = 0
        
        for BlastResult in self.shortBlastResults():
            if BlastResult.Gene == gene:
                hitCounter += 1
            
        return hitCounter
    
    def countBlastHitsByGroup(self, gene, group):
        hitCounter = 0
        
        for BlastResult in self.shortBlastResults():
            if (BlastResult.Gene == gene and BlastResult.NomenclatureFields[0] == group) :
                hitCounter += 1
            
        return hitCounter
        
    def averageScoreByGene(self, gene):
        scoreCounter = 0
        hitCounter = 0
        
        for BlastResult in self.shortBlastResults():
            if BlastResult.Gene == gene:
                hitCounter += 1
                scoreCounter += BlastResult.BlastScore
            
        if hitCounter == 0:
            return 0
        else:
            return (1.0 * scoreCounter) / hitCounter
    
    def averageScoreByGroup(self, gene, group):
        scoreCounter = 0
        hitCounter = 0
        
        for BlastResult in self.shortBlastResults():
            if (BlastResult.Gene == gene and BlastResult.NomenclatureFields[0] == group) :
                hitCounter += 1
                scoreCounter += BlastResult.BlastScore
            
        if hitCounter == 0:
            return 0
        else:
            return (1.0 * scoreCounter) / hitCounter
    
    


