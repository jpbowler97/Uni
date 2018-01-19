import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.backends.backend_pdf import PdfPages


# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')


def analytic_solution_1(nterms, X, T): # Question 2 situation 1 solution
	# A = np.array(nterms)
	W = 0
	for k in range(nterms):
		# print(k)
		# A[k] = 16.0/(np.pow(np.pi, 3)*np.pow(2*k + 1, 3))
		A = 16.0/(math.pow(np.pi, 3)*math.pow(2*k + 1, 3))
		W += A*math.sin((math.pi/2.0 + k*math.pi)*X)*math.exp(-1*(math.pow(math.pi, 2)/4.0 + k*math.pow(math.pi, 2) + math.pow(k*math.pi, 2))*T)
	V = (1.0/2.0)*math.pow(X, 2) - X + W
	U = T + V
	return U

# analytic_solution_1(2, 1, 1)

def analytic_solns_1(nterms, X_values, T_values):
	n = len(X_values)*len(T_values)
	results_table = np.ndarray(shape=(n, 3))
	i = 0
	for t in T_values:
		for x in X_values:
			results_table[i] = (x, t, analytic_solution_1(nterms, x, t))
			i += 1
	# print(results_table)
	df = pd.DataFrame(results_table)
	df.rename(columns={0: 'X', 1: 'T', 2: 'U'}, inplace=True)
	# print(df.head())
	df = df.to_latex()
	afile = open(r'/Users/James/Documents/catam/LaTeX/2.1-TheDiffusionEquation/Q2_0_Results.tex', "w")
	afile.write(df)

	return results_table.transpose()[2] # U values


def analytic_solution_2(nterms, X, T): # Question 2, situation 2 solution
	W = 0
	for k in range(1, nterms+1):
		# print(k)
		B = (-2.0/(3.0*math.pi*k))*(math.pow(-1.0, k)) + (4.0 + 6.0*math.pow(-1.0, k))*(1.0/(math.pow(k*math.pi, 3)))
		W += B*math.sin((k*math.pi)*X)*math.exp(-1*(math.pow(k*math.pi, 2))*T)
	V = (1.0/2.0)*math.pow(X, 2) - (1.0/6.0)*math.pow(X, 3) - (1.0/3.0)*X + W
	U = V + (1 - X)*T
	return U	

def analytic_solns_2(nterms, X_values, T_values):
	n = len(X_values)*len(T_values)
	results_table = np.ndarray(shape=(n, 3))
	i = 0
	for t in T_values:
		for x in X_values:
			results_table[i] = (x, t, analytic_solution_2(nterms, x, t))
			i += 1
	# print(results_table)
	df = pd.DataFrame(results_table)
	df.rename(columns={0: 'X', 1: 'T', 2: 'U'}, inplace=True)
	# print(df.head())
	df = df.to_latex()
	afile = open(r'/Users/James/Documents/catam/LaTeX/2.1-TheDiffusionEquation/Q2_1_Results.tex', "w")
	afile.write(df)
	# print(results_table.transpose())
	return results_table.transpose()[2] # U values

# analytic_solns_2(3, (1, 2, 3), (1, 2, 3))

def analytic_solution_3(X, T): # Question 1 solution
	xi = X/math.sqrt(T)
	F = (1 + (1.0/2.0)*(xi**2))*math.erfc((1.0/2.0)*xi) - math.pow(math.pi, -(1.0/2.0))*xi*math.exp(-((xi**2)/4.0))
	U = T*F
	return U

def analytic_solns_3(X_values, T_values):
	n = len(X_values)*len(T_values)
	results_table = np.ndarray(shape=(n, 3))
	i = 0
	for t in T_values:
		for x in X_values:
			results_table[i] = (x, t, analytic_solution_3(x, t))
			i += 1
	# print(results_table)
	df = pd.DataFrame(results_table)
	df.rename(columns={0: 'X', 1: 'T', 2: 'U'}, inplace=True)
	# print(df.head())
	df = df.to_latex()
	afile = open(r'/Users/James/Documents/catam/LaTeX/2.1-TheDiffusionEquation/Q2_2_Results.tex', "w")
	afile.write(df)

	return results_table.transpose()[2] # U values


def heatflux1(T):
	return (-T - (1.0/2.0)*math.pow(T, 1.0/2.0) + math.pow(math.pi, -1.0/2.0)*math.pow(T, 1.0/2.0))

def heatflux2(T):
	return (1 - (8.0/math.pow(math.pi, 2))*math.exp((-math.pow(math.pi, 2)/4.0)*T))

def heatflux3(T):
	return (T + 1.0/3.0)

	# 3D plot??
	# x = X_values
	# t = T_values
	# X, T = np.meshgrid(x, t)
	# zs = np.array([analytic_solution_1(nterms, x,t) for x,t in zip(np.ravel(X), np.ravel(T))])
	# U = zs.reshape(X.shape)

	# ax.plot_surface(X, T, U)

	# ax.set_xlabel('T Label')
	# ax.set_ylabel('X Label')
	# ax.set_zlabel('U Label')

	# plt.show()

# analytic_solns_1(3, (1, 2, 3), (1, 2))

def Q2():
	nterms = 5
	# T_values = np.array([0.08, 0.24, 0.48, 0.96, 1.92])
	T_values = np.array([0.08])
	X_values = np.arange(0.1, 1.1, 0.1)
	# print(X_values)
	U1 = analytic_solns_1(nterms, X_values, T_values) # also outputs tables
	U2 = analytic_solns_2(nterms, X_values, T_values)
	U3 = analytic_solns_3(X_values, T_values)

	# Issues with U2???

	# plt.plot(X_values, U1, U2, U3)
	with PdfPages('/Users/James/Documents/catam/LaTeX/2.1-TheDiffusionEquation/Q2_NonDim1_8.pdf') as pdf:
		plt.plot(X_values, U3)
		plt.xlabel('X')
		plt.ylabel('U')
		pdf.savefig()
		# plt.show()
		plt.close()

	with PdfPages('/Users/James/Documents/catam/LaTeX/2.1-TheDiffusionEquation/Q2_NonDim2_8.pdf') as pdf:
		plt.plot(X_values, U1)
		plt.xlabel('X')
		plt.ylabel('U')
		pdf.savefig()
		# plt.show()
		plt.close()

	with PdfPages('/Users/James/Documents/catam/LaTeX/2.1-TheDiffusionEquation/Q2_NonDim3_8.pdf') as pdf:
		plt.plot(X_values, U2)
		plt.xlabel('X')
		plt.ylabel('U')
		pdf.savefig()
		# plt.show()
		plt.close()

	T_values = np.array([0.24])
	X_values = np.arange(0.1, 1.1, 0.1)
	# print(X_values)
	U1 = analytic_solns_1(nterms, X_values, T_values) # also outputs tables
	U2 = analytic_solns_2(nterms, X_values, T_values)
	U3 = analytic_solns_3(X_values, T_values)

	# Issues with U2???

	# plt.plot(X_values, U1, U2, U3)
	with PdfPages('/Users/James/Documents/catam/LaTeX/2.1-TheDiffusionEquation/Q2_NonDim1_24.pdf') as pdf:
		plt.plot(X_values, U3)
		plt.xlabel('X')
		plt.ylabel('U')
		pdf.savefig()
		# plt.show()
		plt.close()

	with PdfPages('/Users/James/Documents/catam/LaTeX/2.1-TheDiffusionEquation/Q2_NonDim2_24.pdf') as pdf:
		plt.plot(X_values, U1)
		plt.xlabel('X')
		plt.ylabel('U')
		pdf.savefig()
		# plt.show()
		plt.close()

	with PdfPages('/Users/James/Documents/catam/LaTeX/2.1-TheDiffusionEquation/Q2_NonDim3_24.pdf') as pdf:
		plt.plot(X_values, U2)
		plt.xlabel('X')
		plt.ylabel('U')
		pdf.savefig()
		# plt.show()
		plt.close()

	T_values = np.array([0.48])
	X_values = np.arange(0.1, 1.1, 0.1)
	# print(X_values)
	U1 = analytic_solns_1(nterms, X_values, T_values) # also outputs tables
	U2 = analytic_solns_2(nterms, X_values, T_values)
	U3 = analytic_solns_3(X_values, T_values)

	# Issues with U2???

	# plt.plot(X_values, U1, U2, U3)
	with PdfPages('/Users/James/Documents/catam/LaTeX/2.1-TheDiffusionEquation/Q2_NonDim1_48.pdf') as pdf:
		plt.plot(X_values, U3)
		plt.xlabel('X')
		plt.ylabel('U')
		pdf.savefig()
		# plt.show()
		plt.close()

	with PdfPages('/Users/James/Documents/catam/LaTeX/2.1-TheDiffusionEquation/Q2_NonDim2_48.pdf') as pdf:
		plt.plot(X_values, U1)
		plt.xlabel('X')
		plt.ylabel('U')
		pdf.savefig()
		# plt.show()
		plt.close()

	with PdfPages('/Users/James/Documents/catam/LaTeX/2.1-TheDiffusionEquation/Q2_NonDim3_48.pdf') as pdf:
		plt.plot(X_values, U2)
		plt.xlabel('X')
		plt.ylabel('U')
		pdf.savefig()
		# plt.show()
		plt.close()

	T_values = np.array([0.96])
	X_values = np.arange(0.1, 1.1, 0.1)
	# print(X_values)
	U1 = analytic_solns_1(nterms, X_values, T_values) # also outputs tables
	U2 = analytic_solns_2(nterms, X_values, T_values)
	U3 = analytic_solns_3(X_values, T_values)

	# Issues with U2???

	# plt.plot(X_values, U1, U2, U3)
	with PdfPages('/Users/James/Documents/catam/LaTeX/2.1-TheDiffusionEquation/Q2_NonDim1_96.pdf') as pdf:
		plt.plot(X_values, U3)
		plt.xlabel('X')
		plt.ylabel('U')
		pdf.savefig()
		# plt.show()
		plt.close()

	with PdfPages('/Users/James/Documents/catam/LaTeX/2.1-TheDiffusionEquation/Q2_NonDim2_96.pdf') as pdf:
		plt.plot(X_values, U1)
		plt.xlabel('X')
		plt.ylabel('U')
		pdf.savefig()
		# plt.show()
		plt.close()

	with PdfPages('/Users/James/Documents/catam/LaTeX/2.1-TheDiffusionEquation/Q2_NonDim3_96.pdf') as pdf:
		plt.plot(X_values, U2)
		plt.xlabel('X')
		plt.ylabel('U')
		pdf.savefig()
		# plt.show()
		plt.close()

	T_values = np.array([1.92])
	X_values = np.arange(0.1, 1.1, 0.1)
	# print(X_values)
	U1 = analytic_solns_1(nterms, X_values, T_values) # also outputs tables
	U2 = analytic_solns_2(nterms, X_values, T_values)
	U3 = analytic_solns_3(X_values, T_values)

	# plt.plot(X_values, U1, U2, U3)
	with PdfPages('/Users/James/Documents/catam/LaTeX/2.1-TheDiffusionEquation/Q2_NonDim1_192.pdf') as pdf:
		plt.plot(X_values, U3)
		plt.xlabel('X')
		plt.ylabel('U')
		pdf.savefig()
		# plt.show()
		plt.close()

	with PdfPages('/Users/James/Documents/catam/LaTeX/2.1-TheDiffusionEquation/Q2_NonDim2_192.pdf') as pdf:
		plt.plot(X_values, U1)
		plt.xlabel('X')
		plt.ylabel('U')
		pdf.savefig()
		# plt.show()
		plt.close()

	with PdfPages('/Users/James/Documents/catam/LaTeX/2.1-TheDiffusionEquation/Q2_NonDim3_192.pdf') as pdf:
		plt.plot(X_values, U2)
		plt.xlabel('X')
		plt.ylabel('U')
		pdf.savefig()
		# plt.show()
		plt.close()

	# plots of heat flux
	T_values = np.array([0.08, 0.24, 0.48, 0.96, 1.92])
	hf1 = np.vectorize(heatflux1)
	hf2 = np.vectorize(heatflux2)
	hf3 = np.vectorize(heatflux3)


	with PdfPages('/Users/James/Documents/catam/LaTeX/2.1-TheDiffusionEquation/Q2_HeatFlux_1.pdf') as pdf:
		plt.plot(T_values, hf1(T_values))
		plt.xlabel('T')
		plt.ylabel('Heat Flux')
		pdf.savefig()
		plt.show()
		plt.close()

	with PdfPages('/Users/James/Documents/catam/LaTeX/2.1-TheDiffusionEquation/Q2_HeatFlux_2.pdf') as pdf:
		plt.plot(T_values, hf2(T_values))
		plt.xlabel('T')
		plt.ylabel('Heat Flux')
		pdf.savefig()
		plt.show()
		plt.close()

	with PdfPages('/Users/James/Documents/catam/LaTeX/2.1-TheDiffusionEquation/Q2_HeatFlux_3.pdf') as pdf:
		plt.plot(T_values, hf3(T_values))
		plt.xlabel('T')
		plt.ylabel('Heat Flux')
		pdf.savefig()
		plt.show()
		plt.close()

Q2()

def Numerical_Scheme(C, N, T_upper):
	delta_X = 1/(1.0*N)
	print('delta_X = ' + str(delta_X))
	delta_T = C*(delta_X**2)
	print('delta_T = ' + str(delta_T))
	M = math.ceil(T_upper/delta_T)
	print('max time = ' + str(M*delta_T))
	U = np.zeros(shape=(N+1, M+1)) # U(n, m) = U(n delta_X, m delta_T)
	for m in range(M+1):
		U[0, m] = m*delta_T
		# print('U[0, m] = ' + str(U[0]))
	# use arange?
	# print(U)
	for n in range(N+1): # unnecessary as already initialised as zeros
		U[n, 0] = 0
	# U[0] = np.zeros(M)
	# print(U)
	for m in range(1, M):
		for n in range(1, N):
			if n == N: #(Insulated end conditions (18))
				U[N, m] = U[N-2, m]
			U[n, m+1] = U[n, m] + C*(U[n+1, m] - 2*U[n, m] + U[n-1, m])
	time_1 = int((0.24/delta_T) + 1) # 0.24
	# print(time_1)
	# print(str(U[0, time_1]))
	t1_profile = U.transpose()[time_1]

	time_2 = int((0.24/delta_T) + 1) # 0.24
	time_3 = int((0.48/delta_T) + 1) # 0.48
	time_4 = int((0.96/delta_T) + 1) # 0.96
	time_5 = int((1.92/delta_T) + 1) # 1.92

	t2_profile = U.transpose()[time_2]
	t3_profile = U.transpose()[time_3]
	t4_profile = U.transpose()[time_4]
	t5_profile = U.transpose()[time_5]

	nterms = 10
	X_values = np.arange(0, 1.0 + delta_X, delta_X)
	print('X_values = ' + str(X_values))

	# analytic_solns_1(nterms, X_values, T_values):
	t1_analytic_profile = analytic_solns_1(nterms, X_values, [0.0])
	t2_analytic_profile = analytic_solns_1(nterms, X_values, [0.24])
	t3_analytic_profile = analytic_solns_1(nterms, X_values, [0.48])
	t4_analytic_profile = analytic_solns_1(nterms, X_values, [0.96])
	t5_analytic_profile = analytic_solns_1(nterms, X_values, [1.92])

	T_values = np.array([0.0, 0.24, 0.48, 0.96, 1.92])
	print('times = ' + str(T_values))

	results_table = np.zeros(shape=(len(T_values)*len(X_values), 5))
	
	# T_values
	temp_arr = np.array([])
	for i in range(len(T_values)):
		temp_arr = np.hstack((temp_arr, np.full(len(X_values), T_values[i]))) # T_values
	results_table.transpose()[0] = temp_arr
	# print(temp_arr)

	# X_values
	temp_arr = np.tile(X_values, len(T_values))
	# print(temp_arr)
	results_table.transpose()[1] = temp_arr

	# Analytic_values
	temp_arr = np.hstack((t1_analytic_profile, t2_analytic_profile, t3_analytic_profile, t4_analytic_profile, t5_analytic_profile))
	results_table.transpose()[2] = temp_arr

	# Numerical_values
	temp_arr = np.hstack((t1_profile, t2_profile, t3_profile, t4_profile, t5_profile))
	results_table.transpose()[3] = temp_arr

	# Difference
	results_table.transpose()[4] = abs(results_table.transpose()[3] - results_table.transpose()[2])

	# print(results_table)
	df = (pd.DataFrame(results_table)
			.rename(columns={0: 'T', 1: 'X', 2: 'Analytic Solns', 3: 'Numerical Solns', 4: 'Error'}))
	print(df.head())

	df = df.to_latex()
	afile = open(r'/Users/James/Documents/catam/LaTeX/2.1-TheDiffusionEquation/Q3_0_Results.tex', "w")
	afile.write(df)

	# print(U)

# Numerical_Scheme(0.5, 5, 1.92)


