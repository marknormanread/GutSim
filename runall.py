import os

for i in range(250):
	os.system("python2.7 launch.py -m " + str(i) + " -d results/24-3wkNewGut -dc 24-24 -end 504.0 -seed " + str(i) + " -f -s 30 -lei -v")