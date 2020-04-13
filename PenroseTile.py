import math
import random
import numpy as np
import matplotlib.pyplot as plt

import string
from PIL import Image

from matplotlib.backends.backend_pdf import PdfPages
# matplotlib.use('Agg')


import matplotlib.path as mpath
import matplotlib.patches as mpatches

from Point import Point


class PenroseTile:
	
	
	# Constructor
	def __init__(self, pathName, multiGridDim, k, shiftConstant, scalingFactor, customShift, tileType, shiftVector):
		self.pltPoints = True
		self.pltCrcPrj = False
		
		plt.clf()
		self.pathName = pathName
		
		
		self.multiGridDim = multiGridDim
		self.scalingFactor = scalingFactor
		self.shiftVector = [0 for n in range(self.multiGridDim)]
		# shiftConstant represents the constraint:
		# sum of all shift values = shiftConstant
		# shiftConstant = 1 for "true" Penrose Tiling
		self.shiftConstant = shiftConstant
		self.tileType = tileType
		# Shift vector (gamma_0, gamma_1, ... , gamma_multiGridDim)

		self.shiftVector = shiftVector
		self.genPathPngName()
		# self.normalizeShiftVector()
	
		
		if(customShift):
			self.shiftVector = [1.5561347099679837,
								0.1221553217085172,
								0.6295885218365656,
								0.6038134174473296,
								1.093308029039604]
		
		
		# genNormVectors scales norm vectors by scalingFactor
		self.genNormVectors()
		
		# k is the number of positive parallel hyperplanes on each grid
		self.k = k
		# n is the total number of parallel hyperplanes on each grid
		# n = k pos. hyperplanes + k neg. hyperplanes + 1 plane
		# this 1 plane is the 0th plane that is shift
		self.numHyperplanes = 2*self.k + 1
		
		# numHyperpoints is the total number of other hyperplanes that intersect each hyperplane
		self.numHyperpoints = (self.multiGridDim-1)*self.numHyperplanes
		self.genMultiGrid()
		# Constant to create null graphs
		self.spaceConst = 2

		
	def genMultiGrid(self):
		p = self.getPoint(0, 1, 0, 0)
		self.multiGrid = [[[[p for b in range(-self.k, self.k+1)] for s in range(self.multiGridDim)] for a in range(-self.k, self.k+1)] for r in range(self.multiGridDim)]
	
	def genHyperpoints(self, pltLines):
		xVals = []
		yVals = []
		plt.subplot(1, 2, 1)


		self.fig = plt.figure()
		ax1 = self.fig.add_subplot(1, 2, 1)

		
		# Plot gridlines
		if (pltLines):
			x = np.arange(-(self.k)-2, self.k+2, 0.1)
			for i in range(self.multiGridDim):
				for t in range(-self.k, self.k+1):
					p_i_Kt_hpEqn = self.getHyperplaneEqn(i, t)
					# p_rs_ab
					eqn = x*p_i_Kt_hpEqn[0] + p_i_Kt_hpEqn[1]
					plt.plot(eqn, x, color = 'r', linewidth = 0.5)
					ax1.plot(eqn, x, 'r')
		
		
		# Traverse through sets of two grids
		for r in range(self.multiGridDim):
			for a in range(-self.k, self.k+1):
				for s in range(r+1, self.multiGridDim):
					for b in range(-self.k, self.k+1):
						# p_rs_ab
						pointVect = self.getPoint(r, s, a, b)
						p_rs_a_b = Point(r, s, a, b, pointVect, self.shiftVector[r], self.shiftVector[s])
						# GENERATE TILE VERTICES
						self.genTileVertices(p_rs_a_b)
						self.multiGrid[r][a][s][b] = p_rs_a_b
						
						# Add to lists from plotting
						xVals.append(pointVect[0])
						yVals.append(pointVect[1])
						
						printAllPoints = False
						if(printAllPoints):
							# Print equation of g_r_a
							print("r: " + str(r) + ", a: " + str(a))
							print(str(self.multiGrid[r][a][s][b].r_eqn[0]) + "*x + "  + str(self.multiGrid[r][a][s][b].r_eqn[1]))
							# Print equation of g_s_b
							print("s: " + str(s) + ", b: " + str(b))
							print(str(self.multiGrid[r][a][s][b].s_eqn[0]) + "*x + " + str(self.multiGrid[r][a][s][b].s_eqn[1]))
							print("*"*50)
		
		# Plot hyperpoints
		#################
		if(self.pltPoints):
			plt.scatter(xVals, yVals, s = 1, color = 'blue')
			plt.axis('equal')
			plt.xlim(-self.k - self.spaceConst, self.k + self.spaceConst)
			plt.ylim(-self.k - self.spaceConst, self.k + self.spaceConst)
			plt.title("D: " + str(self.multiGridDim) + ", k: " + str(self.k))
			# plt.show()
		
		# Plot circular projections
		#################
		if(self.pltCrcPrj):
			plt.subplot(1, 2, 1)
			plt.scatter(xVals, yVals, s = 50, color = 'purple', alpha = 0.05)
			plt.axis('equal', square=True)
			plt.xlim(-self.k - self.spaceConst, self.k + self.spaceConst)
			plt.ylim(-self.k - self.spaceConst, self.k + self.spaceConst)
			plt.title("Hyperpoints projected as a function of a constant radius, k: " + str(self.k))
			plt.show()
	
		ax1.scatter(xVals, yVals, s = 1, color = 'blue')
						
	def getPoint(self, r, s, a, b):
		r_ang = r*(math.pi*(2/self.multiGridDim))
		s_ang = s*(math.pi*(2/self.multiGridDim))
		
		sin_r_ang = math.sin(r_ang)
		sin_s_ang = math.sin(s_ang)
		cos_r_ang = math.cos(r_ang)
		cos_s_ang = math.cos(s_ang)
		
		r_a_const = a + 0.5 - self.shiftVector[r]
		s_b_const = b + 0.5 - self.shiftVector[s]
		
		x_rs_ab_num = r_a_const*sin_s_ang - s_b_const*sin_r_ang
		x_rs_ab_den = cos_r_ang*sin_s_ang - cos_s_ang*sin_r_ang
		x_rs_ab = x_rs_ab_num/x_rs_ab_den
		
		y_rs_ab_num = (r_a_const/cos_r_ang) - (s_b_const/cos_s_ang)
		y_rs_ab_den = math.tan(r_ang) - math.tan(s_ang)
		y_rs_ab = y_rs_ab_num/y_rs_ab_den
		
		p_rs_ab = [x_rs_ab, y_rs_ab]
		return p_rs_ab
		
	def getPoint2(self, r, s, a, b):
		x_num = a*self.normVectors[s][1] - b*self.normVectors[r][1]
		x_num+= 0.5*self.normVectors[s][1] - 0.5*self.normVectors[r][1]
		x_num+= self.shiftVector[s]*self.normVectors[r][1] - self.shiftVector[r]*self.normVectors[s][1]
		x_den = self.normVectors[r][0]*self.normVectors[s][1] - self.normVectors[s][0]*self.normVectors[r][1]
		x_rs = x_num / x_den
		
		y_num = a*self.normVectors[s][0] - b*self.normVectors[r][0]
		y_num+= 0.5*self.normVectors[s][0] - 0.5*self.normVectors[r][0]
		y_num+= self.shiftVector[s]*self.normVectors[r][0] - self.shiftVector[r]*self.normVectors[s][0]
		y_den = self.normVectors[r][1]*self.normVectors[s][0] - self.normVectors[s][1]*self.normVectors[r][0]
		y_rs = y_num / y_den
		
		p_rs = (x_rs, y_rs)
		return p_rs
		
	def plotTiles(self, pltLines, showTiles):
		self.ax = plt.subplot(1, 2, 2)
		ax2 = self.fig.add_subplot(1, 2, 2)
		
		# Plot gridlines
		if (pltLines):
			x = np.arange(-(self.k)-2, self.k+2, 0.1)
			for i in range(self.multiGridDim):
				for t in range(-self.k, self.k+1):
					p_i_Kt_hpEqn = self.getHyperplaneEqn(i, t)
					# p_rs_ab
					eqn = x*p_i_Kt_hpEqn[0] + p_i_Kt_hpEqn[1]
					self.ax.plot(eqn, x, color = 'r', linewidth = 0.5)
		
		# Plot tiles
		for r in range(self.multiGridDim):
			for a in range(-self.k, self.k+1):
				for s in range(r+1, self.multiGridDim):
					for b in range(-self.k, self.k+1):
						p_rs_ab = self.multiGrid[r][a][s][b]
						Path = mpath.Path
						# Can be optimized by using for loop to
						# iterate through tile vertices
						path_data = [
							(Path.MOVETO, p_rs_ab.K_p[0]),
							(Path.LINETO, p_rs_ab.K_p[1]),
							(Path.LINETO, p_rs_ab.K_p[2]),
							(Path.LINETO, p_rs_ab.K_p[3]),
							(Path.CLOSEPOLY, p_rs_ab.K_p[0])
						]
						codes, verts = zip(*path_data)
						path = mpath.Path(verts, codes)
						patch = mpatches.PathPatch(path, facecolor='m', lw=1, alpha=0.5)
						self.ax.add_patch(patch)
						
						x, y = zip(*path.vertices)
						# lines, is like that    no syntax error here
						# lines, = self.ax.plot(x, y, 'k')
						ax2.plot(x, y, 'k')

		self.setZero()

		# Configure and show plot		
		self.ax.grid()
		self.ax.axis('equal')
		self.ax.set_xlim(-self.k-2 + self.zero[0], self.k+2 + self.zero[0])
		self.ax.set_ylim(-self.k-2 + self.zero[1], self.k+2 + self.zero[1])
		plt.title("shiftConstant: " + str(self.shiftConstant) + ",  method: " + str(self.tileType))
		
		if(showTiles):
			plt.show()
		
		
		# ax = plt.gca()
		# xVals, yVals = [], []
		# for line in ax.lines:
		# 	x, y = line.get_xdata(), line.get_ydata()
		# 	xVals.append(x)
		# 	yVals.append(y)
		# print(xVals)
		# print(yVals)
		

		# Save image locally
		savingFigure = True
		if(savingFigure):
			self.fig.savefig(self.tilingPathPngName)
		

		plt.close("all")
		
		printTileStats = False
		if(printTileStats):
			print("shift constant: " + str(self.shiftConstant))
			print("gamma: " + self.shiftVectorsToStr())
			print("Swap p_xi&p_yi: " + str(self.swap_x_y) + ", replace K_pi w/ r, s: " + str(self.replace))
			if(self.tileType == "Effinger"):
				print("Effinger method: math.ceil(p_xi + p_yi + self.shiftVector[i])")
			elif(self.tileType == "generlized"):
				print("Generalized method round(p_xi + p_yi - self.shiftVector[i])")
			
		
		
		
	
	
	
	
	
	
	def genTileVertices(self, p):
		K_p = []
		
		for i in range(self.multiGridDim):
			p_xi = (p.p_x)*(self.normVectors[i][0])
			p_yi = (p.p_y)*(self.normVectors[i][1])
			K_pi = 0
			
			self.swap_x_y = False
			if(self.swap_x_y):
				temp = p_xi
				p_xi = p_yi
				p_yi = temp 
			self.replace = False
			if (self.replace):
				if (i == p.r):
					K_pi = p.r
				elif (i == p.s):
					K_pi = p.s
			elif(self.tileType == "Effinger"):
				K_pi = math.ceil(p_xi + p_yi + self.shiftVector[i])
			elif(self.tileType == "roundingSum"):
				K_pi = round(p_xi + p_yi + self.shiftVector[i])
			elif(self.tileType == "generalized"):
				K_pi = round(p_xi + p_yi - self.shiftVector[i])
			elif(self.tileType == "ceilSub"):
				K_pi = math.ceil(p_xi + p_yi - self.shiftVector[i])
			
			K_p.append(K_pi)
		
		x_e_K, y_e_K = 0, 0
		for i in range(self.multiGridDim):
			x_e_K += (self.normVectors[i][0])*(K_p[i])
			y_e_K += (self.normVectors[i][1])*(K_p[i])
		e_K = (x_e_K, y_e_K)
		
		
		
		
		
		
		
		
		x_e_K_r = x_e_K + self.normVectors[p.r][0]
		y_e_K_r = y_e_K + self.normVectors[p.r][1]
		e_K_r = (x_e_K_r, y_e_K_r)
		
		x_e_K_rs = x_e_K_r + self.normVectors[p.s][0]
		y_e_K_rs = y_e_K_r + self.normVectors[p.s][1]
		e_K_rs = (x_e_K_rs, y_e_K_rs)
		
		x_e_K_s = x_e_K + self.normVectors[p.s][0]
		y_e_K_s = y_e_K + self.normVectors[p.s][1]
		e_K_s = (x_e_K_s, y_e_K_s)
		
		p.setTileVertices([e_K, e_K_r, e_K_rs, e_K_s])

	
	def plotHyperplanes(self):
		x = np.arange(-(self.k)-2, self.k+2, 0.1)
		self.p = plt.subplot()
		for i in range(self.multiGridDim):
			for t in range(-self.k, self.k+1):
				p_i_Kt_hpEqn = self.getHyperplaneEqn(i, t)
				# p_rs_ab
				eqn = x*p_i_Kt_hpEqn[0] + p_i_Kt_hpEqn[1]
				#print("i: " + str(i) + " t: " + str(t) + " eqn: " + str(p_i_Kt_hpEqn[0]) + "*x + " + str(p_i_Kt_hpEqn[1]))
				self.p.plot(eqn, x, color = 'r', linewidth=0.5)
			# print("*"*50)
		self.p.set_title("D: " + str(self.multiGridDim) + ", k: " + str(self.k), va='bottom')
		self.p.set_xlim(-self.k - self.spaceConst, self.k + self.spaceConst)
		self.p.set_ylim(-self.k - self.spaceConst, self.k + self.spaceConst)
		plt.show()
				
	def getHyperplaneEqn(self, i, t):
		# Calculate intersection angle
		intAng = i*(math.pi)*(2/self.multiGridDim)
		# Calculate y-intercept
		deltaDiff = t + 0.5 - self.shiftVector[i]
		deltaDiff = deltaDiff / math.cos(intAng)
		# Calculate slope
		arcOperator = -(math.sin(intAng) / math.cos(intAng))
		# hyperplaneEqn = [slope, y-intercept]
		hyperplaneEqn = [arcOperator, deltaDiff]
		return hyperplaneEqn
		
	def genShiftVector(self):
		shiftSum = 0
		# Create random shift values
		for i in range(self.multiGridDim):
			# random() generates shift values over [0.0, 1.0)
			randomInt = random.random()
			self.shiftVector[i] = randomInt
			shiftSum += randomInt
		scaleFactor = shiftSum / self.shiftConstant
		# Re-scale shift values
		for i in range(self.multiGridDim):
			self.shiftVector[i] = self.shiftVector[i] / scaleFactor
			
	def genNormVectors(self):
		self.normVectors = []
		for i in range(self.multiGridDim):
			normVector = self.getNormVector(i)
			self.normVectors.append(normVector)
			
	def getNormVector(self, i):
		normVectorConstant = 2*(math.pi)*(i/self.multiGridDim)
		e_x = math.cos(normVectorConstant)
		e_y = math.sin(normVectorConstant)
		normVector = (self.scalingFactor*e_x, self.scalingFactor*e_y)
		return normVector
		
	def normalizeShiftVector(self):
		shiftSum = 0
		for i in range(self.multiGridDim):
			shiftSum += self.shiftVector[i]
		scaleFactor = shiftSum / self.shiftConstant
		for i in range(self.multiGridDim):
			if scaleFactor == 0:
				self.shiftVector[i] = 0
			else:
				self.shiftVector[i] = self.shiftVector[i]/scaleFactor
			
	def shiftVectorsToStr(self):
		gamma = ""
		gamma = "[ " + str(self.shiftVector[0])
		for i in range(1, self.multiGridDim):
			gamma = gamma + ", "
			gamma = gamma + str(self.shiftVector[i])
		gamma = gamma + "]"
		return gamma
		
	def setShiftVectorOverBound(self, a, b):
		for i in range(self.multiGridDim):
			randomInt = ((a-b)*random.random()) + a
			self.shiftVector[i] = randomInt
		
	def genPathPngName(self):
		self.identifier = "id_" + str(random.randrange(0, 100000)) + "_"
		imageName = str(self.tileType) + '_'
		imageName += 'sC_' + str(self.shiftConstant) + '_'
		imageName += str(self.shiftVector[0])
		for i in range(1, self.multiGridDim):
			imageName += '_'
			shift_i = str(self.shiftVector[i])
			imageName += shift_i
		imageName.replace('.', '')
		self.tilingFigName = imageName
		self.indentifiedFigName = self.identifier + self.tilingFigName
		self.tilingPngName = self.indentifiedFigName + ".png"
		self.tilingPathPngName = self.pathName + self.tilingPngName

	def setNumWhitePixels(self):
		# Open image
		tilingPng = Image.open(self.tilingPathPngName)
		self.tilingImgWidth, self.tilingImgHeight = tilingPng.size

		# Boundaries for second image?
		self.picStart_x = 353
		self.picStart_y = 60

		self.picEnd_x = 574
		self.picEnd_y = 425

		self.numWhitePixels = 0
		for x in range(self.picStart_x, self.picEnd_x +1):
			for y in range(self.picStart_y, self.picEnd_y+1):
				pix = tilingPng.getpixel((x, y))
				pixAve = (pix[0] + pix[1] + pix[2] + pix[3])/4
				if((pixAve >= 250) and (pixAve < 260)):
					self.numWhitePixels += 1

	def setNumWhitePixels2(self, buf):
		numWhitePixels = 0
		for index, value in np.ndenumerate(buf):
			index, value = index, value # Void warn
			numWhitePixels += 1
		

	def setShiftVector(self, sV):
		self.shiftVector = sV
	
	def setZero(self):
		centerPoint = self.multiGrid[0][0][1][0]
		self.zero = centerPoint.K_p[0]


	# def fig2data (self, fig):
	# 	"""
	# 	@brief Convert a Matplotlib figure to a 4D numpy array with RGBA channels and return it
	# 	@param fig a matplotlib figure
	# 	@return a numpy 3D array of RGBA values
	# 	"""
	# 	# draw the renderer
	# 	fig.canvas.draw ( )
	
	# 	# Get the RGB buffer from the figure
	# 	w,h = fig.canvas.get_width_height()


	# 	# For ARGB NParray
	# 	# buf = np.fromstring( fig.canvas.tostring_argb(), dtype=npuint8)
	# 	#buf.shape = (w, h, 4)
	# 	# canvas.tostring_argb give pixmap in ARGB mode. Roll the ALPHA channel to have it in RGBA mode
	# 	# buf = np.roll ( buf, 3, axis = 2 )


	# 	buf = np.fromstring ( fig.canvas.tostring_rgb(), dtype=np.uint8 )
	# 	buf.shape = ( w, h,  3)

	# 	# print(fig.canvas.tostring_rgb())


	# 	# fig.canvas.print_to_buffer()
	
	# def fig2img (self, fig ):
	# 	"""
	# 	@brief Convert a Matplotlib figure to a PIL Image in RGBA format and return it
	# 	@param fig a matplotlib figure
	# 	@return a Python Imaging Library ( PIL ) image
	# 	"""
	# 	# put the figure pixmap into a numpy array
	# 	buf = self.fig2data ( fig )
	# 	w, h, d = buf.shape
	# 	return Image.fromstring( "RGB", ( w ,h ), buf.tostring( ) )