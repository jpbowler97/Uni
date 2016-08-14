# check with other people for accuracy of results to help debug

from math import exp

class ODE:

	# initial_x = 1.0
	# initial_y = 2.0

	def __init__(self, initial_x, initial_y):
		self.initial_x = initial_x
		self.initial_y = initial_y

	def equation(self, x, y): # f = dy/dx = equation(x, y)
		f = -4*y + 3*exp(-x) # for example
		return f

	def Euler_Method(self): # Uses the Euler scheme to evaluate the ODE at a specified x value based on a step_size and nSteps from the initial conditions
		
		# step_size
		step_size = 0.1
		h = step_size

		# nSteps is the number of steps from x_0
		nSteps = 2
		n = nSteps

		# initial conditions
		x_0 = self.initial_x
		y_0 = self.initial_y

		# x-coordinate
		x_n = x_0 + n*h
		print "x-coordinate equals " + str(x_n) 

		# initial x, y values
		x = x_0
		Y = y_0

		# recurssive definitions
		for k in range(1,n): # (i.e. 1 to n-1 inclusive) 
			f = self.equation(x, Y)
			x = x_0 + k*h
			Y += h*f

		return Y

	def Leapfrog_Method(self):

		# step_size
		step_size = 0.1
		h = step_size

		# nSteps is the number of steps from x_0
		nSteps = 2
		n = nSteps

		# initial conditions
		x_0 = self.initial_x
		y_0 = self.initial_y		

		# initial x, y values
		x = x_0
		Y_current = y_0

		# Euler scheme for first value
		f = self.equation(x, Y_current)
		Y_next = Y_current + h*f
		Y_before = Y_current
		Y_current = Y_next

		for k in range(2,n): # nIterations ISSUE: check this works as expected; starts at 2 as first iteration is above 
			f = self.equation(x, Y_current)
			Y_next = Y_before + 2*h*f
			Y_before = Y_current
			Y_current = Y_next
			x = x_0 + k*h

		return Y_current

	def Runge_Kutta_Method(self):

		# step_size
		step_size = 0.1
		h = step_size

		# nSteps is the number of steps from x_0
		nSteps = 2
		n = nSteps

		# initial conditions
		x_0 = self.initial_x
		y_0 = self.initial_y		

		# initial x, y values
		x = x_0
		Y = y_0

		for k in range(1,n):
			f = self.equation(x, Y)
			k_1 = h*f
			f = self.equation(x + (1/2)*h, Y + (1/2)*k_1)
			k_2 = h*f
			f = self.equation(x + (1/2)*h, Y + (1/2)*k_2)
			k_3 = h*f
			f = self.equation(x + h, Y + k_3)
			k_4 = h*f
			Y += (1/6)*(k_1 +2*k_2 + 2*k_3 + k_4)
			x = x_0 + k*h

		return Y

