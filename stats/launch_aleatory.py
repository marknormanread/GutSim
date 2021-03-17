# Will create the directory structure required for, and launch, an aleatory analysis.
#
# Call this script from WITHIN a directory such as exp-aleatory. Otherwise the filenames
# will mess up.
#


import os

def write_pbs(filename, replicates, startSeed, dataPath, mouse, aNum):
	'''
	filename - where the pbs file is to be written.
	replicates - how many replicates, as a STRING.
	startSeed - where seeds should start for this run. As a STRING.
	mouse - which mouse to execute. As a STRING.
	'''
	endReplicates = str(int(replicates) - 1)
	f = open(filename, 'w')
	f.write('#~/bin/sh -login\n')
	f.write('#PBS -l walltime=02:30:00\n')
	f.write('#PBS -q workq\n')
	f.write('#PBS -W group_list=fi3\n')
	f.write('\n')
	f.write('#PBS -N Al' + replicates + '.' + aNum + '\n')
	f.write('#PBS -e OUTPUT.err\n')
	f.write('#PBS -o OUTPUT.out\n')
	f.write('\n')
	f.write('#### array job ####\n')
	f.write('#PBS -J 0-' + endReplicates + '\n')
	f.write('\n')
	f.write('DIR=' + dataPath + '\n')
	f.write('cd ~/micro\n')
	f.write('python experiment.py -m ' + mouse + ' -sm $PBS_ARRAY_INDEX -seed $(($PBS_ARRAY_INDEX+' + startSeed + ')) -d $DIR -dc 24-3 -end 504.0 -f -s 50 -tag ' + replicates + '\n')
	f.close()


ss = [1,5,10,50,100,250]	# sample sizes to investigate
a_replicates = 20			# how many replicates to perform for each sample size, for statistical comparisons.
topPath = 'aleatory'


# create these directories.
for s in ss:
	dir = str(s)
	if not os.path.exists(dir):
		os.makedirs(dir)

# populate directories with replicate experiments.
for s in ss:
	for a in range(a_replicates):
		dir = str(s) + '/' + str(a)
		if not os.path.exists(dir):
			os.makedirs(dir)

# write PBS scripts into replicate experiment directories.
scripts = []			# locations to scripts to be launched using qsub.
startSeed = 0
for s in ss:
	for a in range(a_replicates):
		dataPath = topPath + '/' + str(s) + '/' + str(a)
		fileName = str(s) + '/' + str(a) + '/submit-' + str(s) + '.' + str(a) + '.pbs'
		write_pbs(filename=fileName, replicates=str(s), startSeed=str(startSeed), dataPath=dataPath, mouse=str(47), aNum=str(a))
		startSeed += a_replicates
		scripts.append(fileName)

for f in scripts:
	print 'launching ' + f
	os.system('qsub ' + f)

