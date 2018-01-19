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

x_0 = 0.0
y_0 = 0.0
step_size = 0.4
total_steps = 10

analytical_solns(x_0, step_size, total_steps)
# print(np.multiply(2, (2, 3)))