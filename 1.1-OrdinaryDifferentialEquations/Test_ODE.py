from ODE import *
from Analytical_Solution import analytical_soln
from Stability import stability, growth
	
x_0 = 0.0
y_0 = 0.0
step_size = 0.4
final_step = 25
# nSteps = 0

ode = ODE(x_0, y_0) # initialises ODE object

# em = ode.Euler_Method(step_size, nSteps)	
# lf = ode.Leapfrog_Method(step_size, nSteps)	
# rk4 = ode.Runge_Kutta_Method(step_size, nSteps)

stability_table = stability(x_0, step_size, final_step, ode)
# print stability_table

growth_table = growth(stability_table)
print growth_table