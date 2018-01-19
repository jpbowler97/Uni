# Programs for 2.4-SimulationsofRandomSamplesfromParametricDistributions CATAM Project

import numpy as np 
import math
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd


def unif():
	x = np.random.uniform()
	return x

def exprv(theta, u): # X = F^{-1}(U) has distribution function F
	x = -(1.0/theta)*math.log(1.0-u)
	return x

# def exprv(theta): 
# 	u = unif()
# 	x = exprv(theta, u)
# 	return x

# y = exprv(2.0, unif())
# print(y)
# x = np.random.exponential(1.0/2.0)
# print(x)

def median(theta):
	x = (1/theta)*math.log(2)
	return x


def lik(param, measured_values):
	param = 1.0*param # int to double
	xbar = np.mean(measured_values)
	n = len(measured_values)
	y = math.pow((math.log(2)/(param)), n)*math.pow((1.0/2.0), n*xbar/param)
	return y

def loglik(param, measured_values):
	y = lik(param, measured_values)
	return math.log(y)

def medianMLE(measured_values):
	# apply median where theta is mean MLE
	thetaMLE = 1.0*len(measured_values)/sum(measured_values)
	return median(thetaMLE)

# y = medianMLE([1.0, 2.0, 3.0])
# print(y)

# def loglik(param, theta, n=1):
# 	x=1
# 	measured_values=[]
# 	for i in range(n):
# 		measured_value = exprv(theta)
# 		measured_values.append(measured_value)
# 		x*=lik(param, measured_value)
# 	medianEst = medianMLE(measured_values)
# 	print('\n' + 'Median MLE Estimate: ' + str(medianEst))
# 	medianAct = median(theta)
# 	print('Actual Median: ' + str(medianAct) + '\n')
# 	return math.log(x)

# y = loglik(1.0, 1.0)
# print(y)

# y = loglik([1.0, 1.0], 1.0)
# print(y)

def plotloglik(m_0, m_1, measured_values): #x_0, x_1 bdry x axis values
	# plot the loglikelihood function over a linspace x axis from about 0.5 to 2.0, for arbitrary theta
	if ((m_0<0) | (m_0>m_1)):
		print('line 71 error caught')
		system.exit(1)
	x_axis = np.linspace(m_0, m_1, num=100) # range of median values
	n = len(x_axis)
	y = np.zeros(n)
	for i in range(n):
		y[i] = loglik(x_axis[i], measured_values)
	plt.plot(x_axis, y)
	plt.xlabel('m')
	plt.ylabel('loglikelihood')
	plt.title('loglikelihood function parametrised by median, m')
	# plt.show()
	return (x_axis, y)

# plotloglik(0.1, 2.0, [1.0, 1.0, 1.0], 1.0, 1.0)



def Q2():
	theta = 1.2
	# print('actual median: ' + str(math.log(2)/(1.0*theta)))
	n = 6
	u = np.zeros(n)
	for i in range(n):
		u[i] = unif()
	# print('uniform rvs: ' + str(u))
	e = np.zeros(n)
	for i in range(n):
		e[i] = exprv(theta, u[i])
	# print('exponential rvs: ' + str(e))
	# print('median mle estimate: ' + str(medianMLE(e)))
	(x, y) = plotloglik(0.1, 2.0, e)
	medMLE = medianMLE(e)
	# print(medMLE)
	# print(x)
	# print(y)
	medACT = math.log(2)/(1.0*theta)
	with PdfPages('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q2_loglikelihood.pdf') as pdf:
		plt.plot(x, y)
		plt.xlabel('m')
		plt.ylabel('loglikelihood')
		plt.title('loglikelihood function parametrised by median, m (n=6)')
		pdf.savefig()
		plt.close()
	tex_file = open('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q2_medACT.tex', "w")
	tex_file.write('actual median: ' + str(medACT))
	tex_file.close()
	tex_file = open('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q2_unifs.tex', "w")
	tex_file.write('uniform rvs: ' + str(u))
	tex_file.close()
	tex_file = open('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q2_exps.tex', "w")
	tex_file.write('exponential rvs: ' + str(e))
	tex_file.close()
	tex_file = open('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q2_medMLE.tex', "w")
	tex_file.write('median mle estimate: ' + str(medMLE))
	tex_file.close()
	tex_file = open('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q2_diff.tex', "w")
	tex_file.write('difference between mle and actual median: ' + str(abs(medMLE - medACT)))
	tex_file.close()

# Q2()

def Q3():
	theta = 1.2
	# print('actual median: ' + str(math.log(2)/(1.0*theta)))
	
	n = 25
	u = np.zeros(n)
	for i in range(n):
		u[i] = unif()
	# print('uniform rvs: ' + str(u))
	e = np.zeros(n)
	for i in range(n):
		e[i] = exprv(theta, u[i])
	# print('exponential rvs: ' + str(e))
	# print('median mle estimate: ' + str(medianMLE(e)))
	(x, y) = plotloglik(0.1, 2.0, e)
	medMLE = medianMLE(e)
	medACT = math.log(2)/(1.0*theta)
	with PdfPages('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q3_1_loglikelihood.pdf') as pdf:
		plt.plot(x, y)
		plt.xlabel('m')
		plt.ylabel('loglikelihood')
		plt.title('loglikelihood function parametrised by median, m (n=25)')
		pdf.savefig()
		plt.close()
	tex_file = open('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q3_1_medACT.tex', "w")
	tex_file.write('actual median: ' + str(medACT))
	tex_file.close()
	tex_file = open('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q3_1_unifs.tex', "w")
	tex_file.write('uniform rvs: ' + str(u))
	tex_file.close()
	tex_file = open('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q3_1_exps.tex', "w")
	tex_file.write('exponential rvs: ' + str(e))
	tex_file.close()
	tex_file = open('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q3_1_medMLE.tex', "w")
	tex_file.write('median mle estimate: ' + str(medMLE))
	tex_file.close()
	tex_file = open('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q3_1_diff.tex', "w")
	tex_file.write('difference between mle and actual median: ' + str(abs(medMLE - medACT)))
	tex_file.close()

	n = 50
	u = np.zeros(n)
	for i in range(n):
		u[i] = unif()
	# print('uniform rvs: ' + str(u))
	e = np.zeros(n)
	for i in range(n):
		e[i] = exprv(theta, u[i])
	# print('exponential rvs: ' + str(e))
	# print('median mle estimate: ' + str(medianMLE(e)))
	(x, y) = plotloglik(0.1, 2.0, e)
	medMLE = medianMLE(e)
	medACT = math.log(2)/(1.0*theta)
	with PdfPages('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q3_2_loglikelihood.pdf') as pdf:
		plt.plot(x, y)
		plt.xlabel('m')
		plt.ylabel('loglikelihood')
		plt.title('loglikelihood function parametrised by median, m (n=50)')
		pdf.savefig()
		plt.close()
	tex_file = open('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q3_2_medACT.tex', "w")
	tex_file.write('actual median: ' + str(medACT))
	tex_file.close()
	tex_file = open('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q3_2_unifs.tex', "w")
	tex_file.write('uniform rvs: ' + str(u))
	tex_file.close()
	tex_file = open('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q3_2_exps.tex', "w")
	tex_file.write('exponential rvs: ' + str(e))
	tex_file.close()
	tex_file = open('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q3_2_medMLE.tex', "w")
	tex_file.write('median mle estimate: ' + str(medMLE))
	tex_file.close()
	tex_file = open('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q3_2_diff.tex', "w")
	tex_file.write('difference between mle and actual median: ' + str(abs(medMLE - medACT)))
	tex_file.close()

	n = 100
	u = np.zeros(n)
	for i in range(n):
		u[i] = unif()
	# print('uniform rvs: ' + str(u))
	e = np.zeros(n)
	for i in range(n):
		e[i] = exprv(theta, u[i])
	# print('exponential rvs: ' + str(e))
	# print('median mle estimate: ' + str(medianMLE(e)))
	(x, y) = plotloglik(0.1, 2.0, e)
	medMLE = medianMLE(e)
	medACT = math.log(2)/(1.0*theta)
	with PdfPages('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q3_3_loglikelihood.pdf') as pdf:
		plt.plot(x, y)
		plt.xlabel('m')
		plt.ylabel('loglikelihood')
		plt.title('loglikelihood function parametrised by median, m (n=100)')
		pdf.savefig()
		plt.close()
	tex_file = open('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q3_3_medACT.tex', "w")
	tex_file.write('actual median: ' + str(medACT))
	tex_file.close()
	tex_file = open('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q3_3_unifs.tex', "w")
	tex_file.write('uniform rvs: ' + str(u))
	tex_file.close()
	tex_file = open('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q3_3_exps.tex', "w")
	tex_file.write('exponential rvs: ' + str(e))
	tex_file.close()
	tex_file = open('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q3_3_medMLE.tex', "w")
	tex_file.write('median mle estimate: ' + str(medMLE))
	tex_file.close()
	tex_file = open('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q3_3_diff.tex', "w")
	tex_file.write('difference between mle and actual median: ' + str(abs(medMLE - medACT)))
	tex_file.close()

# Q3()


def dist2rv(theta): #see LaTeX doc for description
	u = unif()
	# print('u: ' + str(u))
	y = (1-u)/math.exp(1)
	# print('y: ' + str(y))
	# we seek the value of k s.t. ke^{-k} = g
	# we use the interval bisection algorithm, with the assumption k < 20
	ub = 20.0
	lb = 1.0
	def g(k):
		return k*math.exp(-k)
	for iters in range(15):
		# print('ub: ' + str(ub))
		# print('lb: ' + str(lb))
		k = (ub + lb)/2.0
		# print('k: ' + str(k))
		# print('g(k): ' + str(g(k)))
		if g(k) < y:
			ub = k
		else:
			lb = k
	# k = theta*x + 1
	x = (k - 1)/theta
	return x

def loglik2(theta, measured_values):
	# see LaTeX for derivation of likelihood function
	s = sum(measured_values)
	p = np.prod(measured_values)
	n = len(measured_values)
	# print(theta)
	# print(p)

	y = 2*n*math.log(theta) + math.log(p) - theta*s
	# print(y)
	return y

def plotloglik2(theta_0, theta_1, measured_values): #x_0, x_1 bdry x axis values
	# plot the loglikelihood function over a linspace x axis from about 0.5 to 2.0, for arbitrary theta
	if ((theta_0<=0) | (theta_0>theta_1)):
		print('line 290 error caught')
		system.exit(1)
	x_axis = np.linspace(theta_0, theta_1, num=100) # range of median values
	n = len(x_axis)
	y = np.zeros(n)
	for i in range(n):
		y[i] = loglik2(x_axis[i], measured_values)
	plt.plot(x_axis, y)
	plt.xlabel('theta')
	plt.ylabel('loglikelihood')
	plt.title('loglikelihood function parametrised by theta, m')
	# plt.show()
	return (x_axis, y)

def thetaMLEest(measured_values):
	return 2.0/(np.mean(measured_values))

def Q7():
	theta = 2.2

	n = 10
	measured_values = np.zeros(n)
	for i in range(n):
		measured_values[i] = dist2rv(theta)
	(x, y) = plotloglik2(0.1, 10, measured_values)
	mle_est = thetaMLEest(measured_values)
	theta_act = 2.2
	with PdfPages('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q7_1_loglikelihood.pdf') as pdf:
		plt.plot(x, y)
		plt.xlabel('theta')
		plt.ylabel('loglikelihood')
		plt.title('loglikelihood function parametrised by theta (n=10)')
		pdf.savefig()
		plt.close()

	tex_file = open('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q7_1_mle.tex', "w")
	tex_file.write('theta mle: ' + str(mle_est))
	tex_file.close()
	tex_file = open('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q7_1_rvs.tex', "w")
	tex_file.write('rvs: ' + str(measured_values))
	tex_file.close()
	tex_file = open('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q7_1_diff.tex', "w")
	tex_file.write('difference between mle and actual theta value: ' + str(abs(mle_est - theta_act)))
	tex_file.close()

	n = 30
	measured_values = np.zeros(n)
	for i in range(n):
		measured_values[i] = dist2rv(theta)
	(x, y) = plotloglik2(0.1, 10, measured_values)
	mle_est = thetaMLEest(measured_values)
	theta_act = 2.2
	with PdfPages('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q7_2_loglikelihood.pdf') as pdf:
		plt.plot(x, y)
		plt.xlabel('theta')
		plt.ylabel('loglikelihood')
		plt.title('loglikelihood function parametrised by theta (n=30)')
		pdf.savefig()
		plt.close()

	tex_file = open('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q7_2_mle.tex', "w")
	tex_file.write('theta mle: ' + str(mle_est))
	tex_file.close()
	tex_file = open('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q7_2_rvs.tex', "w")
	tex_file.write('rvs: ' + str(measured_values))
	tex_file.close()
	tex_file = open('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q7_2_diff.tex', "w")
	tex_file.write('difference between mle and actual theta value: ' + str(abs(mle_est - theta_act)))
	tex_file.close()

	n = 50
	measured_values = np.zeros(n)
	for i in range(n):
		measured_values[i] = dist2rv(theta)
	(x, y) = plotloglik2(0.1, 10, measured_values)
	mle_est = thetaMLEest(measured_values)
	theta_act = 2.2
	with PdfPages('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q7_3_loglikelihood.pdf') as pdf:
		plt.plot(x, y)
		plt.xlabel('theta')
		plt.ylabel('loglikelihood')
		plt.title('loglikelihood function parametrised by theta (n=50)')
		pdf.savefig()
		plt.close()

	tex_file = open('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q7_3_mle.tex', "w")
	tex_file.write('theta mle: ' + str(mle_est))
	tex_file.close()
	tex_file = open('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q7_3_rvs.tex', "w")
	tex_file.write('rvs: ' + str(measured_values))
	tex_file.close()
	tex_file = open('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q7_3_diff.tex', "w")
	tex_file.write('difference between mle and actual theta value: ' + str(abs(mle_est - theta_act)))
	tex_file.close()

# Q7()

# print(dist2rv(2.2))

def dist2check():
	theta = 2.2
	y_hist = np.zeros(1000000)
	for i in range(1000000):
		y_hist[i] = dist2rv(theta)
	# plt.hist(y, bins=1000)
	# plt.show()
	x_0 = min(y_hist)
	# print('x_0: ' + str(x_0))
	x_1 = max(y_hist)
	# print('x_1: ' + str(x_1))
	x_axis = np.linspace(x_0, x_1, 1000)
	# print('x_axis: ' + str(x_axis))
	y = np.zeros(1000)
	i=0
	for x in x_axis:
		y[i] = math.pow(theta, 2)*x*math.exp(-1*theta*x)
		i+=1
	# print('y :' + str(y))

	with PdfPages('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q5_rvdist.pdf') as pdf:
		my_hist = plt.hist(y_hist, bins=500, normed=True)
		distfunc = plt.plot(x_axis, y)
		pdf.savefig()
		# plt.show()


# dist2check()

def Q5():
	dist2check()

# Q5()

# x = np.zeros((2, 3))
# x = np.array([[1, 2, 3], [4, 5, 6]])
# print(x[1])

def thetaMLEs(Nsamples, sample_size, theta):
	N = Nsamples
	n = sample_size
	x = np.zeros((N, n)) 	# intialise zero array with dimensions [N x n] (rows are different samples)
	thetas = np.zeros(N)
	for i in range(N):
		for j in range(n):
			x[i, j] = dist2rv(theta)
		sampleMLE = thetaMLEest(x[i])
		thetas[i] = sampleMLE
	return thetas

def plothist(thetas):
	# plot a histogram of the thetas
	# plt.hist(thetas)
	# plt.xlabel('theta estimate')
	# plt.ylabel('frequency')
	# plt.show()
	with PdfPages('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q8_1_histogram.pdf') as pdf:
		plt.hist(thetas)
		plt.xlabel('theta estimate')
		plt.ylabel('frequency')
		plt.title('Histogram of theta estimates (N=200 samples of size n=10)')
		pdf.savefig()
		plt.close()

def Q8():
	N = 200

	n = 10
	theta = 2.2
	thetaESTs = thetaMLEs(N, n, theta)
	with PdfPages('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q8_1_histogram.pdf') as pdf:
		plt.hist(thetaESTs)
		plt.xlabel('theta estimate')
		plt.ylabel('frequency')
		plt.title('Histogram of theta estimates (N=200 samples of size n=10)')
		pdf.savefig()
		plt.close()

	n = 30
	theta = 2.2
	thetaESTs = thetaMLEs(N, n, theta)
	with PdfPages('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q8_2_histogram.pdf') as pdf:
		plt.hist(thetaESTs)
		plt.xlabel('theta estimate')
		plt.ylabel('frequency')
		plt.title('Histogram of theta estimates (N=200 samples of size n=30)')
		pdf.savefig()
		plt.close()

	n = 50
	theta = 2.2
	thetaESTs = thetaMLEs(N, n, theta)
	with PdfPages('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q8_3_histogram.pdf') as pdf:
		plt.hist(thetaESTs)
		plt.xlabel('theta estimate')
		plt.ylabel('frequency')
		plt.title('Histogram of theta estimates (N=200 samples of size n=50)')
		pdf.savefig()
		plt.close()

# Q8()

def normalrv(mu, sigma):
	# generate normal rv
	A = unif()
	B = unif()
	phi = 2*math.pi*A
	V = -2*math.log(1 - B)
	X = mu + sigma*math.sqrt(V)*math.cos(phi)
	return X

def normalcheck():
	x = np.zeros(1000000)
	for i in range(1000000):
		x[i] = normalrv(0, 1)
	with PdfPages('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q10_1_normalhist.pdf') as pdf:
		plt.hist(x, bins=400, normed=True)
		pdf.savefig()
		plt.close()

# normalcheck()

def intervalcalc(measured_values):
	# calculates boundaries for an 80% confidence interval
	n = len(measured_values)
	lb = np.mean(measured_values) - 1.282/math.sqrt(n)
	ub = np.mean(measured_values) + 1.282/math.sqrt(n)
	return (lb, ub)

def intervalcheck(mu, lb, ub):
	if lb < mu < ub:
		return True
	else:
		return False
	# checks to see if mu is in the interval

def Q11():
	mu = 0
	N = 25
	n = 100
	successes = 0
	fails = 0
	x = np.zeros(n)
	ubs = np.zeros(N)
	lbs = np.zeros(N)
	sample_means = np.zeros(N)
	bools = np.zeros(N)

	for i in range(N):
		for j in range(n):
			x[j] = normalrv(mu, 1)
		(lbs[i], ubs[i]) = intervalcalc(x)
		sample_means[i] = np.mean(x)
		if intervalcheck(mu, lbs[i], ubs[i]) == True:
			bools[i] = 1
	
	from collections import OrderedDict

	d = OrderedDict({'sample mean': sample_means,
		 'lower bound': lbs,	
		 'upper bound': ubs,
		 'mean in interval?': bools})
	df = pd.DataFrame(d)
	fails = len(df[df['mean in interval?'] == 0])
	# print(fails)
	tex_file = open('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q11_1_nfails.tex', "w")
	tex_file.write('number of times the interval didn\'t contain the mean: ' + str(fails))
	tex_file.close()
	# print(df)
	df = df.to_latex()
	tex_file = open(r'/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q11_1_Results.tex', "w")
	tex_file.write(df)
	tex_file.close()

def Q12():
	mu = 4
	N = 1000
	n = 50
	successes = 0
	fails = 0
	x = np.zeros(n)
	ubs = np.zeros(N)
	lbs = np.zeros(N)
	sample_means = np.zeros(N)
	bools = np.zeros(N)

	for i in range(N):
		for j in range(n):
			x[j] = normalrv(mu, 1)
		(lbs[i], ubs[i]) = intervalcalc(x)
		sample_means[i] = np.mean(x)
		if intervalcheck(mu, lbs[i], ubs[i]) == True:
			bools[i] = 1
	
	from collections import OrderedDict

	d = OrderedDict({'sample mean': sample_means,
		 'lower bound': lbs,	
		 'upper bound': ubs,
		 'mean in interval?': bools})
	df = pd.DataFrame(d)
	fails = len(df[df['mean in interval?'] == 0])
	return (1.0*fails)/(1.0*1000)

# print(Q12())

def chi_square(dof):
	n = dof
	x = np.zeros(n)
	for i in range(n):
		x[i] = math.pow(normalrv(0, 1), 2)
	y = sum(x)
	return y

def Q13():
	n = 100

	dof = 1
	x = np.zeros(n)
	for i in range(n):
		x[i] = chi_square(dof)
	with PdfPages('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q13_1_a_hist.pdf') as pdf:
		plt.hist(x)
		pdf.savefig()
		plt.close()

	dof = 5
	x = np.zeros(n)
	for i in range(n):
		x[i] = chi_square(dof)
	with PdfPages('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q13_1_b_hist.pdf') as pdf:
		plt.hist(x)
		pdf.savefig()
		plt.close()

	dof = 40
	x = np.zeros(n)
	for i in range(n):
		x[i] = chi_square(dof)
	with PdfPages('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q13_1_c_hist.pdf') as pdf:
		plt.hist(x)
		pdf.savefig()
		plt.close()

	n = 300

	dof = 1
	x = np.zeros(n)
	for i in range(n):
		x[i] = chi_square(dof)
	with PdfPages('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q13_2_a_hist.pdf') as pdf:
		plt.hist(x, bins=20)
		pdf.savefig()
		plt.close()

	dof = 5
	x = np.zeros(n)
	for i in range(n):
		x[i] = chi_square(dof)
	with PdfPages('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q13_2_b_hist.pdf') as pdf:
		plt.hist(x, bins=20)
		pdf.savefig()
		plt.close()

	dof = 40
	x = np.zeros(n)
	for i in range(n):
		x[i] = chi_square(dof)
	with PdfPages('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q13_2_c_hist.pdf') as pdf:
		plt.hist(x, bins=20)
		pdf.savefig()
		plt.close()

	n = 500

	dof = 1
	x = np.zeros(n)
	for i in range(n):
		x[i] = chi_square(dof)
	with PdfPages('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q13_3_a_hist.pdf') as pdf:
		plt.hist(x, bins=40)
		pdf.savefig()
		plt.close()

	dof = 5
	x = np.zeros(n)
	for i in range(n):
		x[i] = chi_square(dof)
	with PdfPages('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q13_3_b_hist.pdf') as pdf:
		plt.hist(x, bins=40)
		pdf.savefig()
		plt.close()

	dof = 40
	x = np.zeros(n)
	for i in range(n):
		x[i] = chi_square(dof)
	with PdfPages('/Users/James/Documents/catam/LaTeX/2.4-SimulationsofRandomSamplesfromParametricDistributions/Q13_3_c_hist.pdf') as pdf:
		plt.hist(x, bins=40)
		pdf.savefig()
		plt.close()

Q13()


	