'''
Created on 09/09/2014

@author: markread

File reads all files in a directory containing a specified prefix, and compiles the data contained therein into a single
table (of equal dimensions) that represents the median values in each cell.
'''

import glob
import numpy as np
import sys


def compile(directory='.'):


	prefix = directory + '/mouse'
	files = glob.glob(prefix + '*')		# find all files in the directory that start with the prefix.


	# Read the contents of files identified above into a large datastructure.
	data = None
	headerString = None
	lines = []			# for consolidating all the individual files into one.
	for i,f in enumerate(files):
		print 'reading file ' + f

		with open(f, 'r') as f1:
			headerString = f1.readline()
			lines.append(f1.readline())

		d = np.genfromtxt(f, dtype=float, delimiter=',', skip_header=1)
		if data is None:
			dim = [len(files)]
			dim.extend(d.shape)
			data = np.zeros(dim)
		data[i,:] = d

	med = np.median(data,axis=0)	# axis 0 is the innermost.

	out = open(directory + '/median-mouse','w')
	out.write(headerString)
	for i in med:
		out.write(str(i) + ',   ')
	out.close()

	out = open(directory + '/all-mouse','w')
	out.write(headerString)
	for l in lines:
		out.write(l)
	out.close()



# hook for calling script from the command line.
if __name__ == '__main__':
	compile()
