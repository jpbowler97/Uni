#1.1 Code

from math import exp
import numpy as np

# print(np.linspace(0, final_step, int(final_step+1/step_size)))

def analytical_soln(x):
	y = exp(-x) - exp(-4.0*x)
	return y

def analytical_solns(x_0, step_size, total_steps):
	x = x_0
	solns=[]

	# for nSteps in np.linspace(0, final_step, int(final_step/step_size)):
	for count in xrange(0, total_steps):
		# print('xrange', count)
		# print "x-coordinate equals " + str(x)
		y = analytical_soln(x)
		# print "y-coordinate equals " + str(y)	
		# solns.append((x, y))
		solns.append(y)
		x += step_size
	return solns
	
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

def stability(x_0, y_0, step_size, total_steps): # as requested in 2.1 Stability
	stability_table = []
	ode = First_Order_ODE(x_0, y_0)

	for nSteps in range(0,total_steps + 1):
		x = x_0 + step_size*nSteps
		error = ode.Leapfrog_Method(step_size, nSteps) - analytical_soln(x)
		stability_table.append( [x, ode.Leapfrog_Method(step_size, nSteps), analytical_soln(x), error] )# constructs table as described in 2.1 Stability: {x_n, Y_n, y(x_n), E_n = Y_n - y(x_n)}

	return stability_table

def growth(stability_table): # makes table of x and abs(error) from stability_table
	growth_table=[]

	for i in range(0, len(stability_table)):
		growth_table.append( ( stability_table[i][0], abs(stability_table[i][3]) ) )
# plot a graph of x against abs(error) in log-log space and fit to find gradient (gamma/growth rate factor)
# bash out some statistical analysis on gamma? e.g. compute variance/find distribution

	return growth_table

def log_growth(stability_table):
	log_growth_table=[]

	for i in range(1, len(stability_table)):
		log_growth_table.append( ( stability_table[i][0], log(abs(stability_table[i][3])) ) )

	return log_growth_table

def plot_error():
	# g.title('My title')
	g('set style data linespoints')

	x_0 = 0
	y_0 = 0
	step_size = 0.4 # h
	total_steps = int(10.0/step_size)
	stability_table = stability(x_0, y_0, step_size, total_steps)
	plot_points = log_growth(stability_table)
	results = stability_table
	linear_fit = stats.linregress(plot_points) # slope, intercept, r_value, p_value, std_err
	# print(linear_fit)
	# g.plot(growth(stability_table))
	g('set terminal pdf')
	g('set output "/Users/James/Documents/catam/LaTeX/1.1-OrdinaryDifferentialEquations/Q1_0.pdf"')
	g('set xlabel "Independent Variable, x_n"')
	g('set ylabel "log(E_n)"')
	g.plot(log_growth(stability_table))
	df = pd.DataFrame(results)
	df = (df.rename(columns={0:'x_n', 1:'Y_n', 2:'y(x_n)', 3:'E_n = Y_n - y(x_n)'})
			# .set_index('x_n'))
			)
 	df = df.to_latex()
 
	afile=open(r'/Users/James/Documents/catam/LaTeX/1.1-OrdinaryDifferentialEquations/Q1_Results.tex', "wb")
	# afile.write(message)
	afile.write(df)

def plot_gammas():
	x_0 = 0
	y_0 = 0
	step_sizes = [0.01, 0.05, 0.1, 0.2, 0.5, 1.0] # perhaps include more values
	gammas = []
	for step_size in step_sizes:
		total_steps = int(10.0/step_size)
		stability_table = stability(x_0, y_0, step_size, total_steps)
		linear_fit = stats.linregress(log_growth(stability_table)) # slope, intercept, r_value, p_value, std_err
		gammas.append(linear_fit[0])
	d1 = Gnuplot.Data(step_sizes, gammas, with_='lp', title = '\gamma against step size')
	g('set terminal pdf')
	g('set xlabel "Step Size, h"')
	g('set ylabel "Error Growth Rate, \gamma"')
	g('set output "/Users/James/Documents/catam/LaTeX/1.1-OrdinaryDifferentialEquations/Q1_1.pdf"')
	g.plot(d1)

def Q1():
	plot_error()
	plot_gammas()

from First_Order_ODE import *
from Analytical_Solution import *
from math import log
from scipy import stats
import Gnuplot
import numpy as np
import pandas as pd


g = Gnuplot.Gnuplot()

def Euler_RK4_Comparison():
	ode = First_Order_ODE(x_0, y_0)
	Euler_solns = ode.Euler_Method_Solns(step_size, total_steps+1)
	RK4_solns = ode.RK4_Method_Solns(step_size, total_steps+1)
	exact_solns = analytical_solns(x_0, step_size, total_steps+1)
	results = [np.linspace(0, 4, 11), Euler_solns, RK4_solns, exact_solns]
	df = pd.DataFrame(results)
	df = (df.T
			.rename(columns={0:'x', 1:'Euler', 2:'RK4', 3:'Exact'}))
	df = df.to_latex()
 
	afile=open(r'/Users/James/Documents/catam/LaTeX/1.1-OrdinaryDifferentialEquations/Q3_Results.tex', "wb")
	# afile.write(message)
	afile.write(df)

	d1 = Gnuplot.Data(Euler_solns, with_="lp lt rgb 'red' ", title = 'Euler')
	d2 = Gnuplot.Data(RK4_solns, with_="lp lt rgb 'green' ", title = 'RK4')
	d3 = Gnuplot.Data(exact_solns, with_="linespoints lt rgb 'blue' lw 1 pt 5", title = 'Exact')
	g('set title "Iterates For Different Schemes"')
	g('set terminal pdf')
	g('set xlabel "x"')
	g('set ylabel "y"')
	g('set output "/Users/James/Documents/catam/LaTeX/1.1-OrdinaryDifferentialEquations/Q3_0.pdf"')
	g.plot(d1, d2, d3)

def Global_Error_Comparison():
	ode = First_Order_ODE(x_0, y_0)
	upper_bdry = 0.4
	error_table=[]
	x=[]
	y_Euler=[]
	y_LF=[]
	y_RK4=[]
	exact_soln = analytical_soln(upper_bdry)

	for k in range(0, 16):
		# print(k)
		total_steps = pow(2, k)
		step_size = upper_bdry/total_steps
		error_table.append((total_steps, log(step_size), # (total_steps, log(step_size), Euler error, LF error, RK4 error)
				log(abs(ode.Euler_Method(step_size, total_steps) - analytical_soln(upper_bdry))),
				log(abs(ode.Leapfrog_Method(step_size, total_steps) - analytical_soln(upper_bdry))),
				log(abs(ode.Runge_Kutta_Method(step_size, total_steps) - analytical_soln(upper_bdry)))))
		x.append(log(step_size))
		df = pd.DataFrame(error_table)
	 	df = (df
			.rename(columns={0:'Total Steps', 1:'log(Step Size)', 2:'log(Euler Error)', 3:'log(LF Error)', 4:'log(RK4 Error)'}))
	 	df = df.to_latex()

		afile=open(r'/Users/James/Documents/catam/LaTeX/1.1-OrdinaryDifferentialEquations/Q4_Results.tex', "wb")
		# afile.write(message)
		afile.write(df)
		y_Euler.append(log(abs(ode.Euler_Method(step_size, total_steps) - exact_soln)))
		y_LF.append(log(abs(ode.Leapfrog_Method(step_size, total_steps) - exact_soln)))
		y_RK4.append(log(abs(ode.Runge_Kutta_Method(step_size, total_steps) - exact_soln)))
		d1 = Gnuplot.Data(x, y_Euler, with_="lp lt rgb 'red' ", title = 'Euler')
		d2 = Gnuplot.Data(x, y_LF, with_="linespoints lt rgb 'yellow' lw 1 pt 5", title = 'LF')		
		d3 = Gnuplot.Data(x, y_RK4, with_="lp lt rgb 'green' ", title = 'RK4')
		g('set title "log(error) vs log(step size)"')
		g('set xlabel "log(h)"')
		g('set ylabel "log(|E_n|)"')
		g('set terminal pdf')
		g('set output "/Users/James/Documents/catam/LaTeX/1.1-OrdinaryDifferentialEquations/Q4_0.pdf"')
		g.plot(d1, d2, d3)

x_0 = 0
y_0 = 0
step_size = 0.4
total_steps = int(4.0/step_size)
# print(total_steps)
def Q3():
	x_0 = 0
	y_0 = 0
	step_size = 0.4
	total_steps = int(4.0/step_size)
	Euler_RK4_Comparison()

def Q4():
	x_0 = 0
	y_0 = 0
	step_size = 0.4
	total_steps = int(4.0/step_size)
	Global_Error_Comparison()

from Analytical_Solution import *
from math import log, sin, cos, sqrt
from scipy import stats
import Gnuplot
import numpy as np
import pandas as pd

g = Gnuplot.Gnuplot()

class Second_Order_ODE:

	def __init__(self, initial_t, (initial_y_1, initial_y_2), gamma, delta, OMEGA, a, omega):
		self.initial_t = initial_t
		self.initial_y_1 = initial_y_1
		self.initial_y_2 = initial_y_2
		self.gamma = gamma
		self.delta = delta
		self.OMEGA = OMEGA
		self.a = a
		self.omega = omega

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

		for k in range(1,n+1):
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

		for nSteps in range(0, total_steps):
			# print(nSteps)
			y = self.order_2_RK4_Method(step_size, nSteps)
			solns.append( (t, y[0]) )# constructs table as described in 2.1 Stability: {x_n, Y_n, y(x_n), E_n = Y_n - y(x_n)}
			# t = t_0 + step_size*nSteps
			t += step_size # is this printing the first case twice?

		return solns

	def analytical_soln_21(self, t):
		gamma = self.gamma
		omega = self.omega
		# y =(1.0/(pow((1 - pow(omega,2)),2) + pow((omega*gamma),2)))*((omega*exp(-gamma*t/2.0)*((((pow(gamma,2)) - 2*(1 - pow(omega,2)))/(sqrt(4 - pow(gamma,2)) ))*sin((sqrt(4 - pow(gamma,2))/2.0)*t) + gamma*cos((sqrt(4 - pow(gamma,2)))/2.0)*t)) +  (1 - pow(omega,2))*sin(omega*t) - omega*gamma*cos(omega*t)) 
		y = (1.0/7.0)*(exp(-1*t/2.0)*(5*sin(sqrt(3)*t/2.0) + sqrt(3)*cos(sqrt(3)*t/2.0)) - 2*sin(sqrt(3)*t) - sqrt(3)*cos(sqrt(3)*t))
		return y

	def analytical_solns_21(self, step_size, total_steps):
		t_0 = self.initial_t
		t = t_0
		solns=[]

		# for nSteps in np.linspace(0, final_step, int(final_step/step_size)):
		for count in xrange(0, total_steps):
			y = self.analytical_soln_21(t)
			# solns.append((x, y))
			solns.append((t, y))
			t += step_size
		return solns

def Q6(): # need to compute analytical solutions and save plot then import
	initial_t = 0
	t_0 = initial_t
	initial_y_1 = 0
	initial_y_2 = 0
	gamma = 1.0
	delta = 0
	OMEGA = 1.0
	a = 1.0
	omega = sqrt(3)
	step_size = 0.4
	# step_size_values = [0.4, 0.2, 0.1]
	upper_bdry = 10.0
	total_steps = int(upper_bdry/step_size)+1
	ode = Second_Order_ODE(initial_t, (initial_y_1, initial_y_2), gamma, delta, OMEGA, a, omega)
	results_table=[]
	results_table.append(ode.order_2_RK4_Method_Solns(step_size, total_steps))
	results_table.append(ode.analytical_solns_21(step_size, total_steps))
	# print(global_error_table)
	global_error_table=[]
	for i in xrange(0,total_steps):
		global_error_table.append((t_0 + i*step_size,np.subtract(results_table[0], results_table[1])[i][1]))
	# print(global_error_table)
	results_table.append(global_error_table)
	# print(np.shape(results_table))
	# print(results_table[0])
	d1 = Gnuplot.Data(results_table[0], with_="lp lt rgb 'red' ", title = 'Numerical Solutions to 2nd Order ODE')
	d2 = Gnuplot.Data(results_table[1], with_="lp lt rgb 'blue'", title = 'Exact Solutions to 2nd Order ODE')
	g('set terminal pdf')
	g('set output "/Users/James/Documents/catam/LaTeX/1.1-OrdinaryDifferentialEquations/Q6_0.pdf"')
	g('set xlabel "Independent Variable, t"')
	g('set ylabel "Iterates, Y_n"')
	g('set title "Comparison of RK4 to Exact Solutions (h = 0.4)"')
	g.plot(d1, d2)
	df1 = (pd.DataFrame(results_table[0])
			.rename(columns={1: "RK4 Iterates"})) 
	df2 = (pd.DataFrame(results_table[1])
			.rename(columns={1: "Analytical Solutions"}))
	df3 = (pd.DataFrame(results_table[2])
			.rename(columns={1: "Global Error E_n = Y_n - y(x_n)"}))
	# print(df1) # something may be wrong with first 2 entries
	print(df2)
	print(df3)
	df = pd.merge(df1, df3, on=[0])
	print(df)
	df = (pd.merge(df, df2, on=[0])
			.rename(columns={0:'t'}))
	# print(df)

 	df = df.to_latex()
	afile=open(r'/Users/James/Documents/catam/LaTeX/1.1-OrdinaryDifferentialEquations/Q6_0_Results.tex', "wb")
	afile.write(df)

	step_size = 0.2
	# step_size_values = [0.4, 0.2, 0.1]
	upper_bdry = 10.0
	total_steps = int(upper_bdry/step_size)+1
	ode = Second_Order_ODE(initial_t, (initial_y_1, initial_y_2), gamma, delta, OMEGA, a, omega)
	results_table=[]
	results_table.append(ode.order_2_RK4_Method_Solns(step_size, total_steps))
	results_table.append(ode.analytical_solns_21(step_size, total_steps))

	global_error_table=[]
	for i in xrange(0,total_steps):
		global_error_table.append((t_0 + i*step_size,np.subtract(results_table[0], results_table[1])[i][1]))
	# print(global_error_table)
	results_table.append(global_error_table)
	# print(np.shape(results_table))
	# print(results_table[0])
	d1 = Gnuplot.Data(results_table[0], with_="lp lt rgb 'red' ", title = 'Numerical Solutions to 2nd Order ODE')
	d2 = Gnuplot.Data(results_table[1], with_="lp lt rgb 'blue'", title = 'Exact Solutions to 2nd Order ODe')
	g('set terminal pdf')
	g('set output "/Users/James/Documents/catam/LaTeX/1.1-OrdinaryDifferentialEquations/Q6_1.pdf"')
	g('set xlabel "Independent Variable, t"')
	g('set ylabel "Iterates, Y_n"')
	g('set title "Comparison of RK4 to Exact Solutions (h = 0.2)"')
	g.plot(d1, d2)
	df1 = (pd.DataFrame(results_table[0])
			.rename(columns={1: "RK4 Iterates"})) 
	df2 = (pd.DataFrame(results_table[1])
			.rename(columns={1: "Analytical Solutions"}))
	df3 = (pd.DataFrame(results_table[2])
			.rename(columns={1: "Global Error E_n = Y_n - y(x_n)"}))
	# print(df1) # something may be wrong with first 2 entries
	print(df2)
	print(df3)
	df = pd.merge(df1, df3, on=[0])
	print(df)
	df = (pd.merge(df, df2, on=[0])
			.rename(columns={0:'t'}))
	# print(df)

 	df = df.to_latex()
	afile=open(r'/Users/James/Documents/catam/LaTeX/1.1-OrdinaryDifferentialEquations/Q6_1_Results.tex', "wb")
	afile.write(df)

	step_size = 0.1
	# step_size_values = [0.4, 0.2, 0.1]
	upper_bdry = 10.0
	total_steps = int(upper_bdry/step_size)+1
	ode = Second_Order_ODE(initial_t, (initial_y_1, initial_y_2), gamma, delta, OMEGA, a, omega)
	results_table=[]
	results_table.append(ode.order_2_RK4_Method_Solns(step_size, total_steps))
	results_table.append(ode.analytical_solns_21(step_size, total_steps))

	global_error_table=[]
	for i in xrange(0,total_steps):
		global_error_table.append((t_0 + i*step_size,np.subtract(results_table[0], results_table[1])[i][1]))
	# print(global_error_table)
	results_table.append(global_error_table)
	# print(np.shape(results_table))
	# print(results_table[0])
	d1 = Gnuplot.Data(results_table[0], with_="lp lt rgb 'red' ", title = 'Numerical Solutions to 2nd Order ODE')
	d2 = Gnuplot.Data(results_table[1], with_="lp lt rgb 'blue'", title = 'Exact Solutions to 2nd Order ODe')
	g('set terminal pdf')
	g('set output "/Users/James/Documents/catam/LaTeX/1.1-OrdinaryDifferentialEquations/Q6_2.pdf"')
	g('set xlabel "Independent Variable, t"')
	g('set ylabel "Iterates, Y_n"')
	g('set title "Comparison of RK4 to Exact Solutions (h = 0.1)"')
	g.plot(d1, d2)
	df1 = (pd.DataFrame(results_table[0])
			.rename(columns={1: "RK4 Iterates"})) 
	df2 = (pd.DataFrame(results_table[1])
			.rename(columns={1: "Analytical Solutions"}))
	df3 = (pd.DataFrame(results_table[2])
			.rename(columns={1: "Global Error E_n = Y_n - y(x_n)"}))
	# print(df1) # something may be wrong with first 2 entries
	print(df2)
	print(df3)
	df = pd.merge(df1, df3, on=[0])
	print(df)
	df = (pd.merge(df, df2, on=[0])
			.rename(columns={0:'t'}))
	# print(df)

 	df = df.to_latex()
	afile=open(r'/Users/James/Documents/catam/LaTeX/1.1-OrdinaryDifferentialEquations/Q6_2_Results.tex', "wb")
	afile.write(df)

def Q7():
	initial_t = 0
	initial_y_1 = 0
	initial_y_2 = 0
	gamma_values = [0.25, 0.5, 1.0, 1.9]
	delta = 0
	OMEGA = 1.0
	a = 1.0
	omega_values = [1.0, 2.0]
	step_size = 0.4
	upper_bdry = 40.0
	results_table=[]
	for omega in omega_values:
		for gamma in gamma_values:
			ode = Second_Order_ODE(initial_t, (initial_y_1, initial_y_2), gamma, delta, OMEGA, a, omega)
			results_table.append(ode.order_2_RK4_Method_Solns(step_size, int(upper_bdry/step_size)+1))
	g('set terminal pdf')
	g('set output "/Users/James/Documents/catam/LaTeX/1.1-OrdinaryDifferentialEquations/Q7_0.pdf"')
	g('set xlabel "Independent Variable, t"')
	g('set ylabel "Iterates, Y_n"')
	g('set title "omega = 1"')
	d1 = Gnuplot.Data(results_table[0], with_="lp lt rgb 'red' ", title = 'gamma = ' + str(gamma_values[0]))
	d2 = Gnuplot.Data(results_table[1], with_="lp lt rgb 'blue' ", title = 'gamma = ' + str(gamma_values[1]))
	d3 = Gnuplot.Data(results_table[2], with_="lp lt rgb 'green' ", title = 'gamma = ' + str(gamma_values[2]))
	d4 = Gnuplot.Data(results_table[3], with_="lp lt rgb 'yellow' ", title = 'gamma = ' + str(gamma_values[3]))	
	g.plot(d1, d2, d3, d4)
	df1 = (pd.DataFrame(results_table[0])
			.rename(columns={1:'gamma = ' + str(gamma_values[0])}))
	df2 = (pd.DataFrame(results_table[1])
			.rename(columns={1:'gamma = ' + str(gamma_values[1])}))
	df3 = (pd.DataFrame(results_table[2])
			.rename(columns={1:'gamma = ' + str(gamma_values[2])}))
	df4 = (pd.DataFrame(results_table[3])
			.rename(columns={1:'gamma = ' + str(gamma_values[3])}))
	df = pd.merge(df1, df2, on=[0])
	df = pd.merge(df, df3, on=[0])
	df = (pd.merge(df, df4, on=[0])
			.rename(columns={0:'t_n'})
			.drop(df.index[0:15])
			.drop(df.index[26:105]))
	# analytical_solutions = [] # see question
	# df = pd.merge(df, analytical_solutions)
 	df = df.to_latex()
	afile=open(r'/Users/James/Documents/catam/LaTeX/1.1-OrdinaryDifferentialEquations/Q7_0_Results.tex', "wb")
	afile.write(df)
	g('set terminal pdf')
	g('set output "/Users/James/Documents/catam/LaTeX/1.1-OrdinaryDifferentialEquations/Q7_1.pdf"')
	g('set xlabel "Independent Variable, t"')
	g('set ylabel "Iterates, Y_n"')
	g('set title "omega = 2"')
	d1 = Gnuplot.Data(results_table[4], with_="lp lt rgb 'red' ", title = 'gamma = ' + str(gamma_values[0]))
	d2 = Gnuplot.Data(results_table[5], with_="lp lt rgb 'blue' ", title = 'gamma = ' + str(gamma_values[1]))
	d3 = Gnuplot.Data(results_table[6], with_="lp lt rgb 'green' ", title = 'gamma = ' + str(gamma_values[2]))
	d4 = Gnuplot.Data(results_table[7], with_="lp lt rgb 'yellow' ", title = 'gamma = ' + str(gamma_values[3]))	
	g.plot(d1, d2, d3, d4)

	df1 = (pd.DataFrame(results_table[4])
			.rename(columns={1:'gamma = ' + str(gamma_values[0])}))
	df2 = (pd.DataFrame(results_table[5])
			.rename(columns={1:'gamma = ' + str(gamma_values[1])}))
	df3 = (pd.DataFrame(results_table[6])
			.rename(columns={1:'gamma = ' + str(gamma_values[2])}))
	df4 = (pd.DataFrame(results_table[7])
			.rename(columns={1:'gamma = ' + str(gamma_values[3])}))
	df = pd.merge(df1, df2, on=[0])
	df = pd.merge(df, df3, on=[0])
	df = (pd.merge(df, df4, on=[0])
			.rename(columns={0:'t_n'})
			.drop(df.index[0:15])
			.drop(df.index[26:105]))

 	df = df.to_latex()
	afile=open(r'/Users/James/Documents/catam/LaTeX/1.1-OrdinaryDifferentialEquations/Q7_1_Results.tex', "wb")
	afile.write(df)

def Q8():
	initial_t = 0
	initial_y_1 = 0
	initial_y_2 = 0
	gamma = 0
	delta_values = [0.25, 0.5, 1.0]
	# delta_values = [0.25]
	OMEGA = 1.0
	a = 1.0
	omega = 1.0
	step_size = 0.25
	upper_bdry = 40.0
	results_table=[]
	for delta in delta_values:
		ode = Second_Order_ODE(initial_t, (initial_y_1, initial_y_2), gamma, delta, OMEGA, a, omega)
		results_table.append(ode.order_2_RK4_Method_Solns(step_size, int(upper_bdry/step_size)+1))
	g('set terminal pdf')
	g('set output "/Users/James/Documents/catam/LaTeX/1.1-OrdinaryDifferentialEquations/Q8_0.pdf"')
	g('set xlabel "Independent Variable, x"')
	g('set ylabel "Iterates, Y_n"')
	d1 = Gnuplot.Data(results_table[0], with_="lp lt rgb 'red' ", title = 'delta = ' + str(delta_values[0]))
	d2 = Gnuplot.Data(results_table[1], with_="lp lt rgb 'blue' ", title = 'delta = ' + str(delta_values[1]))
	d3 = Gnuplot.Data(results_table[2], with_="lp lt rgb 'green' ", title = 'delta = ' + str(delta_values[2]))
	# d4 = Gnuplot.Data(results_table[3], with_="lp lt rgb 'yellow' ", title = 'delta = ' + str(delta_values[3]))	
	# g.plot(results_table)
	g.plot(d1, d2, d3)
	df1 = (pd.DataFrame(results_table[0])
			.rename(columns={1:'delta = ' + str(delta_values[0])}))
	df2 = (pd.DataFrame(results_table[1])
			.rename(columns={1:'delta = ' + str(delta_values[1])}))
	df3 = (pd.DataFrame(results_table[2])
			.rename(columns={1:'delta = ' + str(delta_values[2])}))
	# df4 = (pd.DataFrame(results_table[3])
	# 		.rename(columns={1:'delta = ' + str(delta_values[3])}))
	df = pd.merge(df1, df2, on=[0])
	df = (pd.merge(df, df3, on=[0])
	# df = (pd.merge(df, df4, on=[0])
			.rename(columns={0:'x_n'})
			.drop(df.index[0:15])
			.drop(df.index[26:165]))
	# df = (df.rename(columns={0:'x_n', 1:'delta = ' + str(delta_values[0]), 2:'delta = ' + str(delta_values[1]), 3:'delta = ' + str(delta_values[2])})			)
 	# print(df.head())
 	df = df.to_latex()
	afile=open(r'/Users/James/Documents/catam/LaTeX/1.1-OrdinaryDifferentialEquations/Q8_0_Results.tex', "wb")
	afile.write(df)

	delta = 20
	results_table=[]
	ode = Second_Order_ODE(initial_t, (initial_y_1, initial_y_2), gamma, delta, OMEGA, a, omega)
	results_table.append(ode.order_2_RK4_Method_Solns(step_size, int(upper_bdry/step_size)+1))
	d4 = Gnuplot.Data(results_table[0], with_="lp lt rgb 'yellow' ", title = 'delta = ' + str(delta))	
	g('set terminal pdf')
	g('set output "/Users/James/Documents/catam/LaTeX/1.1-OrdinaryDifferentialEquations/Q8_1.pdf"')
	g('set xlabel "Independent Variable, x"')
	g('set ylabel "Iterates, Y_n"')
	g.plot(d4)
	# print(results_table)
	df = (pd.DataFrame(results_table[0])
			.rename(columns={1:'delta = ' + str(delta)})
			.rename(columns={0:'x_n'}))

 	df = df.head(n=10).to_latex()
	afile=open(r'/Users/James/Documents/catam/LaTeX/1.1-OrdinaryDifferentialEquations/Q8_1_Results.tex', "wb")
	afile.write(df)
	
