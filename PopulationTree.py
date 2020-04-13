from Population import Population

from copy import copy
import os

import random
import decimal

import csv


class PopulationTree:
	

    # Constructor
    def __init__(self, pickupFromStop, rootName, homeCSVPath, fieldNames, popIndex, lastPopGenes, searchResolution, avgFitStopCrit, stdDevStopCrit,
                    multiGridDim, k, scalingFactor, numGenLimit, popSize, mutationRate,
                    poolSize, shiftConstantRange, shiftRange, methods):
        
        # Initializing passed in parameters
        #################################################################################
        self.firstStep = True
        self.firstPickupStep = False

        self.pickupFromStop = pickupFromStop
        if self.pickupFromStop:
            self.firstPickupStep = True

        self.rootName = rootName
        self.homeCSVPath = homeCSVPath

        self.fieldNames = fieldNames

        # searchResolution = 1000 and shiftRange = 10   ==>   ... 9.998, 9.999, 10.000
        self.searchResolution = searchResolution

        # Stopping Criterion
        self.avgFitStopCrit = avgFitStopCrit
        self.stdDevStopCrit = stdDevStopCrit

        # Multigrid Variables
        self.multiGridDim = multiGridDim
        self.k = k
        self.scalingFactor = scalingFactor

        # Population Variables
        self.numGenLimit = numGenLimit
        # self.popSize MUST BE ODD
        self.popSize = popSize
        self.mutationRate = mutationRate
        self.poolSize = poolSize

        # Dualization Variables
        self.shiftConstantRange = shiftConstantRange
        self.shiftRange = shiftRange
        self.methods = methods
        #################################################################################


        # Generate initial population
        self.pltPoints = True
        # Initialize Population Tree containing a history of all populations
        self.popTree = []
        # Instantiate initial population
        self.popIndex = 0

        # Initial population without pickingUp
        if not pickupFromStop:
            rootPathName = self.makeAndReturnDir(self.popIndex)
            self.rootPop = Population(self.pickupFromStop, self.fieldNames, rootPathName, self.popIndex, self.popSize,
                                    self.multiGridDim, self.k, self.scalingFactor,
                                    self.searchResolution, self.methods, self.shiftConstantRange, self.shiftRange)
            # Create penrose tile objects for each population element
            # THIS ALSO EVALUATES FITNESS FUNCTION FOR EACH OBJECT
            self.rootPop.genPopulation()
            # Calculate population fitness
            self.rootPop.calculatePopulationFitness()
            self.popTree.append(self.rootPop)
        # Initial population by last saved population genes
        else:
            self.popIndex = popIndex
            for n in range(self.popIndex):
                n = n # Void warn
                self.popTree.append(None)
            rootPathName = self.makeAndReturnDir(self.popIndex)
            self.rootPop = Population(self.pickupFromStop, self.fieldNames, rootPathName, self.popIndex, self.popSize,
                                    self.multiGridDim, self.k, self.scalingFactor,
                                    self.searchResolution, self.methods, self.shiftConstantRange, self.shiftRange)
            # Create penrose tile objects for each population element
            # THIS ALSO EVALUATES FITNESS FUNCTION FOR EACH OBJECT
            self.rootPop.genPopFromGenes(lastPopGenes)
            # Calculate population fitness
            self.rootPop.calculatePopulationFitness()
            self.popTree.append(self.rootPop)


    def popIter(self):
        # While populatioin has not met stopping creterion
        while (self.popIndex < self.numGenLimit) and (self.popMeetsStoppingCriterion(self.popTree[self.popIndex])):
            self.writePop()
            # Selection and Crossover
            nextPopData = self.selection(self.popTree[self.popIndex])

            # Mutation
            mutatedNextPopData = self.mutation(nextPopData)

            # Update fitness
            pathName = self.makeAndReturnDir(self.popIndex+1)

            nextPop = Population(self.pickupFromStop, self.fieldNames, pathName, self.popIndex + 1, self.popSize, self.multiGridDim, self.k, self.scalingFactor,
                                self.searchResolution, self.methods, self.shiftConstantRange, self.shiftRange)
            
            # Create penrose tile objects for each population element
            # THIS ALSO EVALUATES FITNESS FUNCTION FOR EACH OBJECT
            nextPop.genPopFromGenes(mutatedNextPopData)

            # Calculate population fitness
            nextPop.calculatePopulationFitness()
            self.popTree.append(nextPop)
            # Clear last population in popTree
            self.popTree[self.popIndex] = None

            self.popIndex += 1


    def popMeetsStoppingCriterion(self, pop):
        for p in pop.pop:
            if p.numWhitePixels == 0:
                p.plotTiles(True, True)
                print("found: " + str(p.shiftVector.tostring()))
                print("method: " + str(p.method))
                return False
        return True
            


    def selection(self, pop):
        outputPopData = []


        # Create Pool
        ########################################################################
        pool = []
        poolGenes = []
        poolCount = 0
        while poolCount < self.poolSize:
            r = random.randrange(0, int(pop.totalFitness))
            itSum = r

            # Find which tile should be added to pool
            i = 0
            while itSum > 0:
                itSum = itSum - pop.popFitness[i]
                i += 1
            pool.append(pop.pop[i-1])
            poolGenes.append(pop.popGenes[i-1])
            poolCount += 1
        
        tempPopGenes = poolGenes
        tempPop = pool
        ########################################################################
        
        # Iterate through population and pop off tiles in pairs of two
        while len(outputPopData) < self.popSize:
            t_1_index = len(tempPop)-1
            tempPop.pop()
            t_2_index = len(tempPop)-1
            tempPop.pop()

            t_1_gene = tempPopGenes[t_1_index]
            t_2_gene = tempPopGenes[t_2_index]

            output_t_1_gene = [0, 0, 0, 0, 0, 0, 0]
            output_t_2_gene = [0, 0, 0, 0, 0, 0, 0]

            rand = random.randrange(0, 2)
            if rand ==0:
                output_t_1_gene = t_1_gene
                output_t_2_gene = t_2_gene
            else:
                output_t_1_gene = t_2_gene
                output_t_2_gene = t_1_gene
            
            outputPopData.append(output_t_1_gene)
            outputPopData.append(output_t_2_gene)
        
        return outputPopData
            
    
    def mutation(self, popData):
        for data in popData:
            for i in range(self.multiGridDim + 2):
                if i != 1:
                    # data = data +- data*mutationRate
                    mut = float(float(random.randrange(0, self.mutationRate*self.searchResolution))/self.searchResolution)
                    mutation = ((-1)**random.randrange(0, 2))*data[i]*mut
                    data[i] =  data[i] + mutation
        return popData


    def makeAndReturnDir(self, i):
        subFldrName = "pop" + str(i) + "/"
        pathName = self.rootName + subFldrName
        try:
            os.mkdir(pathName)
        except:
            print(pathName + " exists, dir not made")
        return pathName


    def writePop(self):
        with open(self.homeCSVPath, 'a', newline='') as w:
            thewriter = csv.DictWriter(w, fieldnames=self.fieldNames)
            if not self.pickupFromStop:
                if self.firstStep:
                    thewriter.writeheader()
                    self.firstStep = False
            thewriter.writerow({'popName': self.popTree[self.popIndex].popName, 'popFldrName': self.popTree[self.popIndex].popPath,
                                'popIndex': self.popIndex, 'multiGridDim': self.multiGridDim, 'k': self.k, 'popSize': len(self.popTree[self.popIndex].pop),
                                'popFitness': self.popTree[self.popIndex].popFitness, 'popGenes': self.popTree[self.popIndex].popGenes})

