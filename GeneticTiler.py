from PopulationTree import PopulationTree

import csv
import shutil
import os

import ast


def main():
    pickupFromStop = True

    rootName = "Pop_Tree_Figs/"
    if not pickupFromStop:
        try:
            shutil.rmtree(rootName[:-1])
            os.mkdir(rootName)
        except:
            print("no such directory to be removed")

    homeCSVName = "popTree.csv"
    homeCSVPath = rootName + homeCSVName

    fieldNames = ['popName', 'popFldrName', 'popIndex', 'multiGridDim', 'k', 'popSize', 'popFitness', 'popGenes']

    step = 0.001
    searchResolution = 1/step

    # Stopping Criterion
    avgFitStopCrit = 0
    stdDevStopCrit = 0

    # Multigrid Variables
    multiGridDim = 5
    k = 3
    scalingFactor = 1

    # Population Variables
    numGenLimit = float('inf')
    # self.popSize MUST BE ODD
    popSize = 1024
    mutationRate = 0.05
    poolSize = 10000

    # Dualization Variables
    shiftConstantRange = 5
    shiftRange = 5
    methods = ['generalized', 'Effinger']
    
    currentPopInd = 0
    lastPopGenes = []
    if pickupFromStop:
        lastPopFldr = ""
        with open(homeCSVPath, 'r') as csv_file:
            csv_reader = csv.DictReader(csv_file)
            for line in csv_reader:
                currentPopInd += 1
                lastPopGenes = line['popGenes']
                lastPopFldr = line['popFldrName']

        # Find missing folder after last folder
        lastPopFldr = lastPopFldr[::-1]
        lastPopFldr = lastPopFldr[1:]
        fldrNum = ""
        # Find folder number iteratively
        while lastPopFldr[0] != 'p':
            fldrNum = lastPopFldr[0] + fldrNum
            lastPopFldr = lastPopFldr[1:]
        fldrNum = str(int(fldrNum) + 1)
        lastPopFldr = lastPopFldr[::-1] + fldrNum
        try:
            shutil.rmtree(lastPopFldr)
        except:
            print("last directory already removed")

        lastPopGenes = ast.literal_eval(lastPopGenes)         

    


    # RUN GENETIC ALGORITHM
    popTree = PopulationTree(pickupFromStop, rootName, homeCSVPath, fieldNames, currentPopInd, lastPopGenes, searchResolution,
                    avgFitStopCrit, stdDevStopCrit, multiGridDim, k, scalingFactor, numGenLimit, popSize, mutationRate,
                    poolSize, shiftConstantRange, shiftRange, methods)
    popTree.popIter()

    


if __name__ == "__main__":
	main()
