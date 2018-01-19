import Gnuplot
import numpy as np

x = np.linspace(0, 10, 100) # http://docs.scipy.org/doc/numpy/reference/generated/numpy.linspace.html
y = x**2
z = 2*y

g = Gnuplot.Gnuplot()

d1 = Gnuplot.Data(x, y, z, with_='lp', title = 'd1')
# d2 = Gnuplot.Data(x, y2, with_='l',title='d2')

g('set grid')
# g('set terminal pdf')
g('set xlabel "hello!"')
# g('set output "/Users/James/Documents/catam/python/1.1-OrdinaryDifferentialEquations/example_plot.pdf"')
g('set key left ')

g.plot(d1)