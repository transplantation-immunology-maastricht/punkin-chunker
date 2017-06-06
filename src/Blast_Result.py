class Blast_Result:
    NomenclatureFields = []
    BlastScore = 0
    Gene = ''
    
    def __init__(self):
        self.NomenclatureFields = []
        self.BlastScore = 0
        self.Gene = ''
    
    def parseNomenclatureLine(self, alleleNameLine):
        # "> HLA-A*02:197 (X2 SIM-I2 X3)"
        #print('Parsing Nomenclature Line:' + alleleNameLine)
        nomenclatureFields = []
        
        alleleNameToken = alleleNameLine.split()[1]
        #print('Allele Name:' + alleleNameToken)
        
        hyphenLocation = alleleNameToken.find('-')
        asteriskLocation = alleleNameToken.find('*')
        
        #print('The Hyphen is at position:' + str(hyphenLocation))        
        #print('The Asterisk is at position:' + str(asteriskLocation))
        
        self.Gene = alleleNameToken[hyphenLocation+1:asteriskLocation]
        #print('Gene Found:' + self.Gene)
        
        nomenclatureFieldString = alleleNameToken[asteriskLocation+1:]        
        #print('Nomenclature String:' + nomenclatureFieldString)
        
        nomenclatureFieldTokens = nomenclatureFieldString.split(':')     
        for nomenclatureFieldToken in nomenclatureFieldTokens:
            nomenclatureFields.append(nomenclatureFieldToken)
        
        self.NomenclatureFields = nomenclatureFields
        
    def parseScoreLine(self, scoreLine):
        # It looks like this:
        # " Score = 1110 bits (601),  Expect = 0.0"
        #print ('Parsing a score line:' + scoreLine)        
        line = ' '.join(scoreLine.split())
        #print ('Line:' + line)
        scoreTokens = line.split(' ')        
        scoreToken = scoreTokens[2]
        #print('Found this score token:' + scoreToken)
        self.BlastScore = float(scoreToken)

