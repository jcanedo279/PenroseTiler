from PenroseTile import PenroseTile


import random
import csv


class Population:
    # Constructor
    def __init__(self, pickupFromStop, fieldNames, pathName, popIndex, popSize, multiGridDim, k, scalingFactor, searchResolution, methods, shiftConstantRange, shiftRange):
        self.pickupFromStop = pickupFromStop
        self.fieldNames = fieldNames
        self.popFieldNames = ['tileName', 'multiGridDim', 'k', 'shiftVector', 'shiftConstant', 'methodIndex', 'fitness', 'genes']


        self.popIndex = popIndex

        self.pathName = pathName
        self.popName = "pop" + str(self.popIndex)
        self.popPath = self.pathName
        self.csvName = self.popName + "saveFile" + ".csv"
        self.csvPath = self.popPath + self.csvName

        self.multiGridDim = multiGridDim
        self.k = k
        self.scalingFactor = scalingFactor

        self.searchResolution = searchResolution
        self.methods = methods
        self.shiftConstantRange = shiftConstantRange
        self.shiftRange = shiftRange

        self.popSize = popSize
        self.pop = []

        self.popGenes = []

        self.tO = 1000000



    def genPopulation(self):
        self.tileInd = 0

        fileMethod = 'w'
        if self.pickupFromStop:
            fileMethod = 'a'
        with open(self.csvPath, fileMethod) as csv_file:
            thewriter = csv.DictWriter(csv_file, fieldnames=self.popFieldNames)
            if self.pickupFromStop:
                with open(self.csvPath, 'r') as csv_file:
                    csv_reader = csv.DictReader(csv_file)
                    for line in csv_reader:
                        line = line # Void warn
                        self.tileInd += 1
                self.tileInd += 1
            with open(self.csvPath, 'w', newline='') as w:
                if(not self.pickupFromStop):
                    w.writeheader()

            customShift = False

            # Create popSize number of PenroseTile objects
            self.popFitness = []
            self.totalFitness = 0
            while self.tileInd < self.popSize:

                # Create random shift constand in range(-sCRange, sCRange) w searchResolution
                s = float(float(random.randrange(-self.shiftConstantRange*self.searchResolution,
                    self.shiftConstantRange*self.searchResolution)) / self.searchResolution)
                # Pick between generalized and Effinger method
                m = self.methods[random.randrange(0, 2)]

                # Create random 5-tupple shifts each in range(-sRange, sRange) w searchResolution
                shiftVector = []
                i = 0
                while i < self.multiGridDim:
                    shiftVector.append(float(float(random.randrange(-self.shiftRange*self.searchResolution,
                                        self.shiftRange*self.searchResolution)) / self.searchResolution))
                    i += 1

                tiley = PenroseTile(self.pathName, self.multiGridDim, self.k, s, self.scalingFactor, customShift, m, shiftVector)
                pltLines = False
                tiley.genHyperpoints(pltLines)
                pltLines = False
                tiley.plotTiles(pltLines, False)
                
                # EVALUATE FITNESS
                tiley.setNumWhitePixels()
                self.popFitness.append(self.tO/tiley.numWhitePixels)
                self.totalFitness += self.tO/tiley.numWhitePixels
                self.pop.append(tiley)

                # Create tile genes non-locally with respect to PenroseTile object
                tileGene = [tiley.shiftConstant, self.methods.index(tiley.tileType), tiley.shiftVector[0], tiley.shiftVector[1],
                            tiley.shiftVector[2], tiley.shiftVector[3], tiley.shiftVector[4]]
                self.popGenes.append(tileGene)
                thewriter.writerow({'tileName': tiley.tilingFigName, 'multiGridDim': tiley.multiGridDim,
                                    'k': tiley.k, 'shiftVector': tiley.shiftVector, 'shiftConstant': tiley.shiftConstant, 
                                    'methodIndex': self.methods.index(tiley.tileType), 'fitness': self.popFitness[self.tileInd-1],
                                    'genes': self.popGenes[self.tileInd]})
                self.tileInd += 1


    def genPopFromGenes(self, popData):
        self.tileInd = 0

        with open(self.csvPath, 'a', newline='') as csv_file:
            thewriter = csv.DictWriter(csv_file, fieldnames=self.popFieldNames)
            if(not self.pickupFromStop):
                thewriter.writeheader()

        customShift = False
        # Create popSize number of PenroseTile objects
        self.popFitness = []
        self.totalFitness = 0
        while self.tileInd < len(popData):
            data = popData[self.tileInd]

            # Create random shift constand in range(-sCRange, sCRange) w searchResolution
            s = data[0]
            # Pick between generalized and Effinger method
            m = int(data[1])

            # Create random 5-tupple shifts each in range(-sRange, sRange) w searchResolution
            shiftVector = []
            for j in range(self.multiGridDim):
                shiftVector.append(data[j+2])

            tiley = PenroseTile(self.pathName, self.multiGridDim, self.k, s, self.scalingFactor, customShift, self.methods[m], shiftVector)
            pltLines = False
            tiley.genHyperpoints(pltLines)
            pltLines = False
            tiley.plotTiles(pltLines, False)
            
            # EVALUATE FITNESS
            tiley.setNumWhitePixels()

            self.pop.append(tiley)
            self.popFitness.append(self.tO/tiley.numWhitePixels)
            self.totalFitness += self.tO/tiley.numWhitePixels

            # Create tile genes non-locally with respect to PenroseTile object
            tileGene = [tiley.shiftConstant, self.methods.index(tiley.tileType), tiley.shiftVector[0], tiley.shiftVector[1],
                        tiley.shiftVector[2], tiley.shiftVector[3], tiley.shiftVector[4]]
            self.popGenes.append(tileGene)

            with open(self.csvPath, 'a') as csv_file:
                tw = csv.DictWriter(csv_file, fieldnames=self.popFieldNames)
                tw.writerow({'tileName': tiley.tilingFigName, 'multiGridDim': tiley.multiGridDim,
                                    'k': tiley.k, 'shiftVector': tiley.shiftVector, 'shiftConstant': tiley.shiftConstant, 
                                    'methodIndex': self.methods.index(tiley.tileType), 'fitness': self.popFitness[self.tileInd-1],
                                    'genes': self.popGenes[self.tileInd-1]})
            self.tileInd += 1


    def calculatePopulationFitness(self):
        self.popFitness = []
        totalFitness = 0
        for tile in self.pop:
            totalFitness += self.tO/tile.numWhitePixels
            self.popFitness.append((self.tO/tile.numWhitePixels)**2)
        # POPULATION MEAN STATISTICS
        self.totalFitness = totalFitness**2
        self.avgFitness = self.totalFitness/self.popSize

        populationVariance = 0
        for tile in self.pop:
            populationVariance += (tile.numWhitePixels - self.avgFitness)**2
        # POPULATION VARIANCE STATISTICS
        self.populationVariance = populationVariance
        self.avgPopulationVariance = self.populationVariance/self.popSize
        
        # Sort normFitnesses, pop, popGenes by fitness least to greatest

        if not self.pickupFromStop:
            zipped = zip(self.popFitness, self.pop, self.popGenes)
            zipped = sorted(zipped, key=lambda tup: tup[0])
            self.popFitness, self.pop, self.popGenes = zip(*zipped)
