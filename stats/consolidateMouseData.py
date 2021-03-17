'''
This script finds all the "mouse-xxx" files in a directory, and compiles the contents of those files into one
larger file.
 example call:
  $> python consolidateMouseData.py <directory to process>
'''#

import glob
import sys 		# pass command line arguments
import re


# for sorting alphanumeric text. Taken from
# http://stackoverflow.com/questions/2669059/how-to-sort-alpha-numeric-set-in-python
# on 22/04/2014.
def sort_alphanumeric( l ):
	""" Sort the given iterable in the way that humans expect."""
	convert = lambda text: int(text) if text.isdigit() else text
	alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
	return sorted(l, key = alphanum_key)


directory = sys.argv[1]	# the first argument is the name of the script to be run.

first = True			# used to write the header line, which should be done once only.

# find files containing mouse data, sort files by index (or number).
files = glob.glob(directory + '/mouse-*')
files = [f for f in files if 'timeseries' not in f]
files = sort_alphanumeric(files)
print('found ' + str(len(files)) + ' individual mouse data files.')

fout = open(directory + '/mice-data.txt','w')
for fileName in files:
	f = open(fileName, 'r')
	lines = f.readlines()
	if first:
		fout.write(lines[0])	# write header line to the file
		first = False
	fout.write(lines[1])		# write second line to the file
	f.close()

fout.close()
