# Module contains useful statistical tools and functions. 
#
#

def Atest(dist1, dist2):
	''' 
	Implementation of Vargha-Delaney's A test, a non-parametric
	effect magnitude test. The distributions may be of different
	sizes in this implementation. 
	'''
	equal = 0.0
	greater = 0.0
	for x in dist1:
		for y in dist2:
			if x == y   : equal += 1
			elif x > y 	: greater += 1
	# multiplication of the number of samples in each distribution
	nm = len(dist1) * len(dist2)	
	return (greater / nm) + ((0.5 * equal) / nm)


def ecdf(dist):
	'''
	Calculates the X and Y coorindates for an empirical cumulative
	distribution plot. It is surprisingly difficult to find a 
	good implementation (easy to install and use), so I've written
	my own.

	'dist' contains the empirically observed values. They do not 
	need to be sorted. 
	
	Returns a tuple, the X coordinates and the 'cdf' which contains
	monotonically values between 0 to 1. Both lists are of the same 
	length, and can be used to plot the graph. 
	'''
	X = sorted(dist)	# X coorinates
	cdf = []			# will put the corresponding CDF values in here. 
	n = len(dist)		
	for xi in X:		# scan through each x-coordinate
		smallEq = 0.0	# count how many values in the distribution are smaller
		for j in dist:
			if xi >= j : smallEq += 1
		prop = smallEq / n # calculate proportion
		cdf.append(prop)
	return (X, cdf)