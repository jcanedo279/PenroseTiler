import math

class Point:
		

	# Constructor takes:
	# i   index of grid
	# k   index of hyperplane
	def __init__(self, r, s, a, b, p, g_r, g_s):
		# grid numbers s.t: 0<=r<s<multiGridDim
		self.r = r
		self.s = s
		
		# a is the hyperplane index of grid r
		# b is the hyperplane index of grid s
		self.a = a
		self.b = b
		
		# p is the point (kept here for efficiency)
		self.p = p
		self.p_x = p[0]
		self.p_y = p[1]
		
		self.r_eqn = self.getHyperplaneEqn(r, a, g_r)
		self.s_eqn = self.getHyperplaneEqn(s, b, g_s)
		
	def setTileVertices(self, K_p):
		self.K_p = K_p
		
		
	def getHyperplaneEqn(self, i, t, g_i):
		# Calculate grid angles
		i_ang = i*(math.pi)*(2/5)
		cos_i_ang = math.cos(i_ang)
		sin_i_ang = math.sin(i_ang)
		# Calculate y-intercept
		interceptConstant = t + 0.5 - g_i
		hyperplaneIntercept = interceptConstant / cos_i_ang
		# Calculate slope
		hyperplaneSlope = -(sin_i_ang / cos_i_ang)
		# hyperplaneEqn = [slope, y-intercept]
		hyperplaneEqn = [hyperplaneSlope, hyperplaneIntercept]
		return hyperplaneEqn
		
		
		
		
		
