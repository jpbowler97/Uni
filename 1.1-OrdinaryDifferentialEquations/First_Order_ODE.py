# check with other people for accuracy of results to help debug
import numpy as np
from math import exp, sin, sqrt
import sys

class First_Order_ODE:

	def __init__(self, initial_x, initial_y):
		self.initial_x = initial_x
		self.initial_y = initial_y

	def equation_1(self, x, y): # f = dy/dx = equation_1(x, y)
		f = -4*y + 3*exp(-x) # for example
		return f

	def Euler_Method(self, step_size, nSteps): # Uses the Euler scheme to evaluate the ODE at a specified x value based on a step_size and nSteps from the initial conditions
		
		# step_size
		h = step_size

		# nSteps is the number of steps from x_0
		n = nSteps

		# initial conditions
		x_0 = self.initial_x
		y_0 = self.initial_y

		# x-coordinate
		x_n = x_0 + n*h
		# print "x-coordinate equals " + str(x_n) 

		# initial x, y values
		x = x_0
		Y = y_0

		# recurssive definitions
		for k in range(1,n): # (i.e. 1 to n-1 inclusive) 
			f = self.equation_1(x, Y)
			x = x_0 + k*h
			Y += h*f

		return Y

	def Euler_Method_Solns(self, step_size, total_steps):
		solns = []
		x_0 = self.initial_x
		x = x_0

		for nSteps in range(1, total_steps + 1):
			y = self.Euler_Method(step_size, nSteps)
			# solns.append( (x, y) )# constructs table as described in 2.1 Stability: {x_n, Y_n, y(x_n), E_n = Y_n - y(x_n)}
			solns.append(y)
			x = x_0 + step_size*nSteps

		return solns

	def Leapfrog_Method(self, step_size, nSteps):

		h = step_size
		x_0 = self.initial_x
		y_0 = self.initial_y		
		x = x_0
		Y_current = y_0

		if nSteps > 0:
		# Euler scheme for first value
			f = self.equation_1(x, Y_current)
			Y_next = Y_current + h*f
			Y_before = Y_current
			Y_current = Y_next
			x = x_0 + h

			for k in range(2,nSteps+1): # nIterations ISSUE: check this works as expected; starts at 2 as first iteration is above 
				f = self.equation_1(x, Y_current)
				Y_next = Y_before + 2*h*f
				Y_before = Y_current
				Y_current = Y_next
				x = x_0 + k*h

		# print(Y_current)

		return Y_current

	def Leapfrog_Method_Solns(self, step_size, total_steps):
		solns = []
		x_0 = self.initial_x

		for nSteps in range(1, total_steps + 1):
			x = x_0 + step_size*nSteps
			y = self.Leapfrog_Method(step_size, nSteps)
			# solns.append( (x, y) )# constructs table as described in 2.1 Stability: {x_n, Y_n, y(x_n), E_n = Y_n - y(x_n)}
			solns.append(y)
		return solns

	def Runge_Kutta_Method(self, step_size, nSteps):

		# step_size
		h = step_size
		# print step_size

		# nSteps is the number of steps from x_0
		n = nSteps

		# initial conditions
		x_0 = self.initial_x
		y_0 = self.initial_y		

		# initial x, y values
		x = x_0
		Y = y_0

		for k in range(1,n+1): # EDIT HERE
			f = self.equation_1(x, Y)
			# print(f)
			k_1 = h*f
			# print k_1
			f = self.equation_1(x + (0.5)*h, Y + (0.5)*k_1)
			# print f
			k_2 = h*f
			f = self.equation_1(x + (0.5)*h, Y + (0.5)*k_2)
			k_3 = h*f
			f = self.equation_1(x + h, Y + k_3)
			k_4 = h*f
			x = x_0 + k*h
			Y += (1.0/6.0)*(k_1 +2*k_2 + 2*k_3 + k_4)
			# print Y

		return Y

	def RK4_Method_Solns(self, step_size, total_steps):
		solns = []
		x_0 = self.initial_x
		x = x_0

		for nSteps in range(1, total_steps+1):
			# print(nSteps)
			# print(nSteps)
			y = self.Runge_Kutta_Method(step_size, nSteps)
			# solns.append( (x, y) )# constructs table as described in 2.1 Stability: {x_n, Y_n, y(x_n), E_n = Y_n - y(x_n)}
			solns.append(y)
			x = x_0 + step_size*nSteps

		return solns		





