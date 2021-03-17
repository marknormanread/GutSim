'''
Created on 16/04/2014

@author: markread

Plots probability curves for how bacteria cells respond to internalised limitting resource, in terms of cellular
division, death, and utilisation of resources.

Probabilities are expressed in terms of cell doing X in an hour, given this level of limiting nutrient.
'''
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import numpy as np
import math
from scipy import stats


if False:
	ln = np.arange(0.0,100.0, 1.0)  # The limiting resources
	# division probabilities
	div = np.exp((ln-100) * (1.0/20.0))
	# death probabilities
	die = np.exp(-ln* (1.0/5.0))
	# rate of internalising alternative nutrient sources
	altHarv = np.exp(-ln * (1.0/20.0)) * 0.33

	plt.clf()
	plt.rcParams.update({'font.size': 22})
# 	p1, = plt.plot(ln,div,'b')
# 	p2, = plt.plot(ln,die,'k')
	p3, = plt.plot(ln,altHarv,'r')
	plt.xlabel('Limiting resource level')
	plt.ylabel('Probability')
	plt.ylim((0.0,1.0))
	plt.legend([p3],['alt src'])
	#plt.show()
	plt.savefig('probability-curves/alt-nut-src.png',dpi=600)

# death rates that incorporate a stress response. When resource becomes limitting the death rate increases, up to
# a point, after which the death rate reduces again (to simulate the aquisition of a stress response phenotype).
if False:
	ln = np.arange(0.0,100.0, 1.0)		# the limitting resources
	d = np.exp(-ln * (1.0/5.0))
	s = np.exp(-ln * (1.0/3.0))
	die = d - s
	plt.clf()
	plt.rcParams.update({'font.size': 22})
	# p1, = plt.plot(ln,d,'b')
	# p2, = plt.plot(ln,s,'k')
	p3, = plt.plot(ln,die,'r')
	plt.xlabel('Limiting resource level')
	plt.ylabel('Probability')
	plt.ylim((0.0,1.0))
	plt.legend([p1,p2,p3],['d','s','die (d-s)'])
	#plt.show()
	plt.savefig('probability-curves/death-growth.png', dpi=600)

# Plot death and growth rates on one graph. 
if False:
	ln = np.arange(0.0,100.0, 1.0)		# the limitting resources
	d = np.exp(-ln * (1.0/5.0))
	s = np.exp(-ln * (1.0/3.0))
	die = d - s
	g_logistic = 1.0 / (1.0 + np.exp((-ln + 90.0) * (1.0/20.0)))
	plt.clf()
	plt.rcParams.update({'font.size': 22})
	p2, = plt.plot(ln, g_logistic,'b')
	p3, = plt.plot(ln,die,'r')
	plt.xlabel('Limiting resource level')
	plt.ylabel('Probability')
	plt.ylim((0.0,1.0))
	plt.legend([p3,p2],['death','growth'])
	#plt.show()
	plt.savefig('death-growth.png', dpi=600)

# experimenting with logistics growth curve.
if False:
	ln = np.arange(0.0,120.0, 1.0)		# the limitting resources
	g_logistic = 1.0 / (1.0 + np.exp((-ln + 90.0) * (1.0/20.0)))
	g_exp = np.exp((ln-100) * (1.0/20.0))

	plt.clf()
	plt.rcParams.update({'font.size': 22})
	p1, = plt.plot(ln,g_logistic,'b')
	p2, = plt.plot(ln,g_exp,'g')
	plt.xlabel('Limiting resource level')
	plt.ylabel('Probability')
	plt.legend([p1,p2],('g-logistic','g-exp'))
	plt.savefig('probability-curves/logistic-growth.png', dpi=600)
	#plt.show()


# experimenting with logistic SI absorption rates, given quantity of consumed nutrients. Decided not to work with
# logistic curves, see notebook for details.
if False:
	qn = np.arange(0,7.5,0.01)
	l_p = 0.8 / (1.0 + np.power(0.2,(-qn+1.61))) + 0.2	# protein absorption
	l_c = 0.8 / (1.0 + np.power(0.2,(-qn+1.95))) + 0.2	# carb absorption
	plt.clf()
	plt.rcParams.update({'font.size': 19})
	pp, = plt.plot(qn,l_p,'c')
	pw, = plt.plot(qn,l_c,'b')
	plt.xlabel('Quantity consumed (grams)')
	plt.ylabel('Proportion')
	plt.legend([pp,pw],('protein','carbohydrate'))
	plt.ylim((0.0,1.0))
	plt.grid()
	plt.savefig('probability-curves/si_absorption_proportion.png', dpi=600)
	#plt.show()

	siConsumedProt = qn * l_p
	siConsumedCarb = qn * l_c
	plt.clf()
	pp, = plt.plot(qn,siConsumedProt,'c')
	pw, = plt.plot(qn,siConsumedCarb,'b')
	plt.xlabel('Quantity consumed (grams)')
	plt.ylabel('Quantity absorbed (grams)')
	plt.legend([pp,pw],('protein','carbohydrate'))
	plt.grid()
	plt.savefig('probability-curves/si_absorbed.png', dpi=600)

	colonEntryProt = qn - siConsumedProt
	colonEntryCarb = qn - siConsumedCarb
	plt.clf()
	plt.rcParams.update({'font.size': 19})
	pp, = plt.plot(qn,colonEntryProt,'c')
	pw, = plt.plot(qn,colonEntryCarb,'b')
	plt.xlabel('Quantity consumed (grams)')
	plt.ylabel('Quantity enters colon (grams)')
	plt.legend([pp,pw],('protein','carbohydrate'))
	plt.grid()
	plt.savefig('probability-curves/colon_supply.png', dpi=600)

# experimenting with rates of small intestine absorption of feed nutrients. This uses increasing exponentially decaying
# curves.
if True:
	qn = np.arange(0, 3.0, 0.01)
	m1 = math.e
	m2 = math.e
	cc = 2.0			# the asymptote to be approached
	cd = 2.0			# the asymptote to be approached
	cp = 2.0			# the asymptote to be approached
	kw = -math.log(1.0 - ((0.34 * 0.70)/cc), m1) / 0.34	  # wheat starch carbohydrate
	kd = -math.log(1.0 - ((0.24 * 0.94)/cd), m1) / 0.24	  # destrinised wheat starch
	kp = -math.log(1.0 - ((0.37 * 0.88)/cp), m1) / 0.37	  # protein
	siConsumeWheat = cc * (1.0 - np.power(m1, -kw * qn))
	siConsumeDex =   cd * (1.0 - np.power(m1, -kd * qn))
	siConsumeProt =  cp * (1.0 - np.power(m1, -kp * qn))   # red
	plt.clf()
	pw, = plt.plot(qn, siConsumeWheat, 'b', linewidth=2.0)
	pd, = plt.plot(qn, siConsumeDex, 'g', linewidth=2.0)
	pp, = plt.plot(qn, siConsumeProt, 'r', linewidth=2.0)
	pyx, = plt.plot(qn, qn, 'k', linewidth=2.0)   # y = x, curve should not breach this, as it implies SI absorbs more than is eaten.
	#pyx = plt.plot(qn, 0.9*qn)
	#plt.title('SI Nutrient Absorption Rates')
	plt.xlabel('Quantity consumed (grams/day)')
	plt.ylabel('Quantity absorbed (grams/day)')
	plt.legend([pw, pd, pp, pyx], ('wheatstarch', 'dextrinised cornstarch', 'protein', 'Y=X'), loc='upper left')
	plt.gca().grid(True, linewidth=2.0)   # turn on grid lines.
	font = {'size'   : 18.0}
	plt.rc('font', **font)
	plt.ylim((0.0, 2.5))
	plt.savefig('si_absorbed.png', dpi=600)
	plt.show()


# plot death and growth rate curves against levels of limiting nutrient.
if False:
	ln = np.arange(0.0,120.0, 1.0)		# the limiting resources
	steepness = 1.0 / 20.0
	g_logistic = 1.0 / (1.0 + np.exp( -steepness * (ln - 90.0)))

	d = np.exp(-ln * (1.0/5.0))
	s = np.exp(-ln * (1.0/3.0))
	die = d - s

	plt.clf()	
	p_g, = plt.plot(ln,g_logistic,'b', linewidth=2.0)
	#p_g, = plt.plot(ln,g_logistic2,'g')
	p_d, = plt.plot(ln,die,'r', linewidth=2.0)
	plt.xlabel('Limiting resource level')
	plt.ylabel('Rate (probability/hour)')
	plt.legend([p_g,p_d],('division','death'),loc='upper left')
	plt.gca().grid(True, linewidth=2.0)   # turn on grid lines.
	font = {'size'   : 18.0}
	plt.rc('font', **font)
	plt.savefig('death_growth_rates.png', dpi=600)
	plt.show()


# Devise a different death rate curve that kills less cells and resutls in more entering the resistant phenotype.
# Model this using a gaussian curve.
if False:
	ln = np.arange(0.01,120.0, 0.01)		# the limiting resources
	# new gaussian death rate curve.
	mean = 2.5
	sigma = 0.3		# standard deviation.
	# Draw the probability density function for this log-normal distribution.
	z = np.exp( -np.power((np.log(ln) - mean),2) / (2*sigma*sigma))
	z *= 1 / (ln * sigma * np.sqrt(2* np.pi))
	# scale values along the y-axis, since gaussian distributions integrate to a value of 1.0 (larger variance = lower
	# peak value)
	z *= 1.0
	plt.plot(ln, z, 'k')

# 	M = float(5.0) # Geometric mean == median
# 	s = float(5.0) # Geometric standard deviation
# 	shape = s # Scipy's shape parameter
# 	scale = M # Scipy's scale parameter
# 	pdf = stats.lognorm.pdf(ln, shape, loc=0, scale=scale)
# 	plt.plot(ln,pdf,'k')


	# previous death rate curve.
	g_logistic = 1.0 / (1.0 + np.exp((-ln + 90.0) * (1.0/20.0)))
	d = np.exp(-ln * (1.0/5.0))
	s = np.exp(-ln * (1.0/3.0))
	die = d - s
	plt.plot(ln, die, 'r')

# 	plt.ylim((0.0,0.5))
# 	plt.xlim((0.0,60.0))
	plt.savefig('probability-curves/death_rate_gaussian.png', dpi=600)
	plt.show()


# devise a model for vairable nutrient uptake rates. As limiting resource drops, so too does capacity to internalise
# nutrient
if False:
	ln = np.arange(0.00,120.0, 0.01)		# the limiting resources
	steepness = 1.0 / 20.0
	logistic = 1.0 / (1.0 + np.exp(-steepness * (ln - 40.0)) )

	plt.clf()
	plt.rcParams.update({'font.size': 22})
	p, = plt.plot(ln,logistic,'b')
	plt.ylim((0.0, 1.0))
	plt.xlabel('Limiting resource level')
	plt.ylabel('Uptake Factor')
	#plt.savefig('probability-curves/uptake_factor.png', dpi=600)
	plt.show()