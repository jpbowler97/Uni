# check with other people for accuracy of results to help debug
import numpy as np
from math import exp, sin, sqrt
import sys

class First_Order_ODE:

	def __init__(x):
		if x.length == 2:
			self.initial_x = x[0]
			self.initial_y = x[1]
		elif x.length == 8:
			self.initial_t = x[0]
			self.initial_y_1 = x[1][0]
			self.initial_y_2 = x[1][1]
			self.gamma = x[2]
			self.delta = x[3]
			self.OMEGA = x[4]
			self.a = x[5]
			self.omega = x[6]
		else: sys.exit("Wrong number of arguments provided to constructor")

	def equation_1(self, x, y): # f = dy/dx = equation_1(x, y)
		f = -4*y + 3*exp(-x) # for example
		return f

	def equation_2(self, t, (y_1, y_2)):
		gamma = self.gamma
		delta = self.delta
		OMEGA = self.OMEGA
		a = self.a
		omega = self.omega

		f_1 = y_2
		f_2 = -gamma*y_2 - pow(delta, 3)*pow(y_1, 2)*y_2 - pow(OMEGA, 2)*y_1 + a*sin(omega*t)
		f = (f_1, f_2)
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

		for nSteps in range(0, total_steps + 1):
			y = self.Euler_Method(step_size, nSteps)
			solns.append( (x, y) )# constructs table as described in 2.1 Stability: {x_n, Y_n, y(x_n), E_n = Y_n - y(x_n)}
			x = x_0 + step_size*nSteps

		return solns

	def Leapfrog_Method(self, step_size, nSteps):

		# step_size
		h = step_size

		# nSteps is the number of steps from x_0
		n = nSteps

		# initial conditions
		x_0 = self.initial_x
		y_0 = self.initial_y		

		# initial x, y values
		x = x_0
		Y_current = y_0

		# Euler scheme for first value
		f = self.equation_1(x, Y_current)
		Y_next = Y_current + h*f
		Y_before = Y_current
		Y_current = Y_next

		for k in range(2,n): # nIterations ISSUE: check this works as expected; starts at 2 as first iteration is above 
			f = self.equation_1(x, Y_current)
			Y_next = Y_before + 2*h*f
			Y_before = Y_current
			Y_current = Y_next
			x = x_0 + k*h

		return Y_current

	def Leapfrog_Method_Solns(self, step_size, total_steps):
		solns = []
		x_0 = self.initial_x

		for nSteps in range(0, total_steps + 1):
			x = x_0 + step_size*nSteps
			y = self.Leapfrog_Method(step_size, nSteps)
			solns.append( (x, y) )# constructs table as described in 2.1 Stability: {x_n, Y_n, y(x_n), E_n = Y_n - y(x_n)}

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

		for k in range(1,n):
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
		print('hello')
		return Y

	def RK4_Method_Solns(self, step_size, total_steps):
		solns = []
		x_0 = self.initial_x
		x = x_0

		for nSteps in range(0, total_steps + 1):
			# print(nSteps)
			y = self.Runge_Kutta_Method(step_size, nSteps)
			solns.append( (x, y) )# constructs table as described in 2.1 Stability: {x_n, Y_n, y(x_n), E_n = Y_n - y(x_n)}
			x = x_0 + step_size*nSteps

		return solns		

	def order_2_RK4_Method(self, step_size, nSteps):

		h = step_size
		n = nSteps

		# x, y_1, y_2, gamma, delta, OMEGA, a, omega
		gamma = self.gamma
		delta = self.delta
		OMEGA = self.OMEGA
		a = self.a
		omega = self.omega

		# initial conditions
		t_0 = self.initial_t
		y_1_0 = self.initial_y_1
		y_2_0 = self.initial_y_2		

		# initial x, y values
		t = t_0
		y_1 = y_1_0
		y_2 = y_2_0
		y = (y_1, y_2)

		for k in range(1,n):
			f = self.equation_2(t, y)
			# print(f)
			k_1 = np.multiply(h, f)
			# print k_1
			f = self.equation_2(t + (0.5)*h, y + np.multiply(0.5, k_1))
			# print f
			k_2 = np.multiply(h, f)
			f = self.equation_2(t + (0.5)*h, y + np.multiply(0.5, k_2))
			k_3 = np.multiply(h, f)
			f = self.equation_2(t + h, y + k_3)
			k_4 = np.multiply(h, f)
			t = t_0 + k*h
			y += np.multiply((1.0/6.0),(k_1 + np.multiply(2, k_2) + np.multiply(2, k_3) + k_4))
			# print Y

		return y

	def order_2_RK4_Method_Solns(self, step_size, total_steps):
		solns = []
		t_0 = self.initial_t
		t = t_0

		for nSteps in range(0, total_steps + 1):
			# print(nSteps)
			y = self.order_2_RK4_Method(step_size, nSteps)
			solns.append( (t, y[0]) )# constructs table as described in 2.1 Stability: {x_n, Y_n, y(x_n), E_n = Y_n - y(x_n)}
			t = t_0 + step_size*nSteps

		return solns





