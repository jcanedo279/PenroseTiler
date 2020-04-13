import csv


from PenroseTile import PenroseTile

def genTiles(pathName):

	fieldNames = ['tileName', 'placeInLoop', 'multiGridDim', 'k', 'shiftVector', 'shiftConstant', 'methodIndex', "numWhitePixels"]

	fileMethod = 'w'
	pickupFromStop = True
	calculateMostFit = False
	if(pickupFromStop):
		fileMethod = 'r'
	with open('tileObjects.csv', fileMethod, newline='') as f:

		numTiles = 0

		methods = ["generalized", "Effinger"]
		tO = 0.001
		r = 10
		res = r
		shiftConstant = 10
		sCInit, aInit, bInit, cInit, dInit, eInit = 0, 0, 0, 0, 0, 0
		if(pickupFromStop):
			lastMInd = 0
			lastPIL = [0 for n in range(5)]
			with open('tileObjects.csv', 'r') as csv_file:
				csv_reader = csv.DictReader(csv_file)
				for line in csv_reader:
					numTiles += 1
					lastMInd = int(line['methodIndex'])
					lastPIL = line['placeInLoop']
			# Convert from string to float
			lastPIL = lastPIL[1:-1]
			lastPIL = lastPIL.split(',')

			sCInit = int(lastPIL[0])
			aInit = int(lastPIL[1])
			bInit = int(lastPIL[2])
			cInit = int(lastPIL[3])
			dInit = int(lastPIL[4])
			eInit = int(lastPIL[5])

			fMethod = True
			fSC, fA, fB, fC, fD, fE = True, True, True, True, True, True
			m = 0
			while m < len(methods):
				if m != lastMInd:
					methods = methods[1:]
					continue
				break
		else:
			sCInit = -shiftConstant
			aInit, bInit, cInit, dInit, eInit = -r, -r, -r, -r, -r
			fMethod = False
			fSC, fA, fB, fC, fD, fE = False, False, False, False, False, False
		
		firstStep = True
		firstIt = True
		mInd = 0
		if(calculateMostFit):
			methods = []
		while mInd < len(methods):
			if(fMethod):
				mInd = lastMInd
				fMethod = False
			sC = -shiftConstant
			while sC < shiftConstant+1:
				if(fSC):
					sC = sCInit
					fSC = False
				a = -r
				while a < r+1:
					if(fA):
						a = aInit
						fA = False
					b = -r
					while b < r+1:
						if(fB):
							b= bInit
							fB = False
						c = -r
						while c < r+1:
							if(fC):
								c = cInit
								fC = False
							d = -r
							while d < r+1:
								if(fD):
									d = dInit
									fD = False
								e = -r
								while e < r+1:
									if(fE):
										e = eInit
										fE = False
									if(pickupFromStop):
										if(firstIt):
											firstIt = False
											# Skip current iteration
											if e!= r:
												e += 1
											else:
												d += 1
									# tO = timesOver, meaning how many times larger than r each shift is
									customShift = False
									if(r == 0):
										shiftVector = [a, b, c, d, e]
									else:
										shiftVector = [tO*a/res, tO*b/res, tO*c/res, tO*d/r, tO*e/res]
									
									tiley = PenroseTile(pathName, 5, 3, tO*sC/shiftConstant, 1, customShift, methods[mInd], shiftVector)
									

									# tiley.plotHyperplanes()
									pltLines = True
									tiley.genHyperpoints(pltLines)
									pltLines = False
									tiley.plotTiles(pltLines)

									tiley.setNumWhitePixels()


									pIL = [sC, a, b, c, d, e]

									with open('tileObjects.csv', 'r'):
										with open('tileObjects.csv', 'a', newline='') as w:
											thewriter = csv.DictWriter(w, fieldnames=fieldNames)
											if(firstStep):
												if(not pickupFromStop):
													thewriter.writeheader()
												firstStep = False
											thewriter.writerow({'tileName': tiley.tilingFigName, 'placeInLoop': pIL, 'multiGridDim': tiley.multiGridDim,
															'k': tiley.k, 'shiftVector': tiley.shiftVector, 'shiftConstant': tiley.shiftConstant, 
															'methodIndex': mInd, 'numWhitePixels': tiley.numWhitePixels})
											numTiles += 1
									e += 1
								d += 1
							c += 1
						b += 1
					a += 1
				sC += 1
			mInd += 1

def mostFit():
	fitnessDict = {}

	with open('tileObjects.csv', 'r') as csv_file:
		csv_reader = csv.DictReader(csv_file)
		numLines = 0
		for line in csv_reader:
			numLines += 1
			key = line['tileName']
			val = line['numWhitePixels']

			fitnessDict.update({key:val})
		print(numLines)

		listofTuples = sorted(fitnessDict.items() ,  key=lambda x: x[1])

		fieldNames = ['tileName', "numWhitePixels"]

		with open('tileObjectsSorted.csv', 'w', newline='') as w:
			thewriter = csv.DictWriter(w, fieldnames=fieldNames)
			thewriter.writeheader()

			for elem in listofTuples:
				thewriter.writerow({'tileName': elem[0], "numWhitePixels": elem[1]})

			



















def main():
	pathName = "figs_cache/"
	methods = ["generalized", "Effinger", "roundingSum", "ceilSub"]
	
	customShift = False
	configuration = False
	shiftation = False
	repitition = False
	narrowing = False

	fitness = True


	if(fitness):
		genTiles(pathName)
		mostFit()

	elif(configuration):
		configurations(pathName, methods, customShift)
	
	elif(shiftation):
		shiftations(pathName, customShift)
	
	elif(repitition):
		repititions(pathName, customShift)
	
	elif(narrowing):
		narrowings(pathName, customShift)
	else:
		### TESTING ###
		tiley = PenroseTile(pathName, 5, 6, 0.000001, 1, customShift, "generalized", [])
		tiley.setShiftVectorOverBound(-1, 1)
		
		# tiley.plotHyperplanes()
		pltLines = True
		tiley.genHyperpoints(pltLines)
		pltLines = False
		tiley.plotTiles(pltLines)
		tiley.setNumWhitePixels()

if __name__ == "__main__":
	main()

def configurations(pathName, methods, customShift):
	scalingFactors = [-2, -1, -0.5, -0.00001, 0.00001, 0.5, 1, 2]
	for method in methods:
		for i in scalingFactors:
			tiley = PenroseTile(pathName, 5, 0, i, 1, customShift, method, [])
			pltLines = True
			tiley.genHyperpoints(pltLines)
			pltLines = False
			tiley.plotTiles(pltLines)

def shiftations(pathName, customShift):
	d = 30
	totalSum = 2*d+1
	for i in range (-d, d+1):
		if (i == 0):
			x = 0.000001
		else:
			x = 5*i/totalSum
		tiley = PenroseTile(pathName, 5, 3, x, 1, customShift, "generalized", [])
	
		pltLines = True
		tiley.genHyperpoints(pltLines)
		pltLines = False
		tiley.plotTiles(pltLines)

def repititions(pathName, customShift):
	for i in range(100):
		tiley = PenroseTile(pathName, 5, 6, i, 1, customShift, "generalized", [])

		pltLines = True
		tiley.genHyperpoints(pltLines)
		pltLines = False
		tiley.plotTiles(pltLines)


def narrowings(pathName, customShift):
	d = 100
	totalSum = 2*d+1
	for i in range(-d, d+1):
		x = (i/(5*totalSum))
		tiley = PenroseTile(pathName, 5, 6, x, 1, customShift, "generalized", [])
		
		pltLines = True
		tiley.genHyperpoints(pltLines)
		pltLines = False
		tiley.plotTiles(pltLines)
