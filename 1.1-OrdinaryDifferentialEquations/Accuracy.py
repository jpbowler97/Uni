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



