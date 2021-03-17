# Checks that all 250 mouse simulation output files are present in a particular directory.
#
#



import os.path
import sys      # provides access to command line arguments.

dir = '.'
if len(sys.argv) > 1:
  dir = sys.argv[1]


missing = False
for i in range(250):
	filename = dir + '/mouse-' + str(i)
	if not os.path.isfile(filename):
		print('missing file number ' + str(i))
		missing = True

if not missing:
	print('all files present in directory ' + dir)

