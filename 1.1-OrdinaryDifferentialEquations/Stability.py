from First_Order_ODE import *
from Analytical_Solution import analytical_soln
from math import log
from scipy import stats
import Gnuplot
import numpy as np
import pandas as pd

g = Gnuplot.Gnuplot()


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

