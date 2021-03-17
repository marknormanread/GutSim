# Generate figure for Andy, relating to the effect of fiber on bacterial landscapes
#
#

import numpy as np
import matplotlib.pyplot as plt

# 0 - mouse_ID
# 1 - carb-intake-KJ/d
# 2 - prot-intake-KJ/d
# 3 - relBcf
# 4 - relBcm
# 5 - relBdf
# 6 - relBdm
# 7 - relBmm
# 8 - absBcf
# 9 - absBcm
# 10 - absBdf
# 11 - absBdm
# 12 - absBmm
# 13 - total
# 14 - nLimBcf
# 15 - nLimBcm
# 16 - nLimBdf
# 17 - nLimBdm
# 18 - nLimBmm

def plotRel(name, hours, relBcf, relBcm, relBdf, relBdm, relBmm):
	plt.clf()
	plt.plot(hours,relBcf, 'm-')
	plt.plot(hours,relBcm, 'b-')
	plt.plot(hours,relBdf, 'r-')
	plt.plot(hours,relBdm, 'c-')
	plt.plot(hours,relBmm, 'g-')
	(ymin,ymax) = plt.ylim()
	plt.fill([0,3,3,0], [0,0,ymax,ymax], 'b', alpha=0.1)
	plt.xlim([0,24])
	plt.ylabel('Relative abundance (%)')
	plt.xlabel('Time after meal (hours)')
	plt.savefig(name, dpi=600)

def plotAbs(name, hours, absBcf, absBcm, absBdf, absBdm, absBmm):
	plt.clf()
	plt.plot(hours,absBcf, 'm-')
	plt.plot(hours,absBcm, 'b-')
	plt.plot(hours,absBdf, 'r-')
	plt.plot(hours,absBdm, 'c-')
	plt.plot(hours,absBmm, 'g-')
	(ymin,ymax) = plt.ylim()
	plt.fill([0,3,3,0], [0,0,ymax,ymax], 'b', alpha=0.1)
	plt.xlim([0,24])
	plt.ylabel('Absolute abundance')
	plt.xlabel('Time after meal (hours)')
	plt.savefig(name, dpi=600)

def plotAvg(mice, data):
	(d1, d2, d3) = data.shape
	# will store the median response data for the indicated mice here.
	# first dimension is the hour, the second is the responses.
	medians = np.zeros([d1,d3])
	for i in range(d1):
		# mice is a list, so this selects the responses for the given mice and hour.
		selected = data[i,mice,:]
		avg = np.mean(selected,axis=0)
		medians[i,:] = avg

	plotAbs('ABS_14-57-29-MED',hours,medians[:,8],medians[:,9],medians[:,10],medians[:,11],medians[:,12])
	plotRel('REL_14-57-29-MED',hours,medians[:,3],medians[:,4],medians[:,5],medians[:,6],medians[:,7])




ctrlName = 'ctrl'
#hours = [0,3,9,12,15,18,21,0]
hours = [0,3,6,12,15,0]

ctrl = np.loadtxt(ctrlName,delimiter=',',skiprows=1)
data = None
for i, h in enumerate(hours):
	d = np.loadtxt(str(h),delimiter=',',skiprows=1)
	# if not already initialised, initialise the size of data.
	if data == None:
		(d1,d2) = d.shape
		data = np.zeros([len(hours),d1,d2])
	data[i,:,:] = d

hours[-1] = 24

dietMice = [24,25,74,75,124,125,174,175,224,225]
plotAvg(dietMice, data)
#mouse = 44
if False:
	for mouse in range(250):
		relBcf = [i[mouse,3] for i in data]
		relBcm = [i[mouse,4] for i in data]
		relBdf = [i[mouse,5] for i in data]
		relBdm = [i[mouse,6] for i in data]
		relBmm = [i[mouse,7] for i in data]

		absBcf = [i[mouse,8] for i in data]
		absBcm = [i[mouse,9] for i in data]
		absBdf = [i[mouse,10] for i in data]
		absBdm = [i[mouse,11] for i in data]
		absBmm = [i[mouse,12] for i in data]

		plotRel('relAb_mouse' + str(mouse), hours, relBcf, relBcm, relBdf, relBdm, relBmm)
		plotAbs('AbAb_mouse' + str(mouse), hours, absBcf, absBcm, absBdf, absBdm, absBmm)


