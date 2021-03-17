__author__ = 'mark'

import math


def phi_normal(x):
	"""
	Provides the integration of a standard normal Gaussian (mu=0, sigma=1) from -Inf to x. Makes use of the error
	function to do so.
	"""
	y = x / math.sqrt(2)
	return 0.5 + (0.5 * math.erf(y))


def phi_gauss(x, mu, sigma):
	"""
	Provides the integration of a Gaussian (with any mean and standard deviation), from -Inf to x. Makes use of the
	error function to do so.

	phi(x, mu, sigma) = phi_normal( (x-mu)/sigma )

	:param x: value to integrate to, starting at -Inf
	:param mu: mean of the Gaussian.
	:param sigma: standard deviation of the Gaussian.
	:return:
	"""
	y = (float(x) - float(mu)) / float(sigma)   # make necessary adjustments for non-normal Gaussian.
	return phi_normal(y)


