'''
Created on 16/04/2014

@author: markread
'''
import random
import matplotlib.pyplot as plt
import numpy as np

plt.clf()
n = np.arange(0.0,100.0, 1.0)

## Either of these next two lines works.
div = np.exp(n-100)
div2 = np.exp((n-100) / 2.0)
div3 = (n / 100) ** 2
#div = np.exp(n) / np.exp(100)

die = np.exp(-n)
die2 = np.exp(-n/20.0)
die3 = (-(n/100)) ** 2
die4 = ((100-n)/100) ** 4

if False:
	plt.plot(n,div,'b')
	plt.plot(n,div2,'r')
	plt.plot(n,div3,'g')
	plt.show()

plt.plot(n,die)
plt.plot(n,die2)
plt.plot(n,die3)
plt.plot(n,die4)

#plt.plot(n,div)
plt.show()

#plt.ylim((0,1.0))
#plt.yscale('log')


