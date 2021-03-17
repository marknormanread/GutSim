# Plots the graphs of diet and nutrition analysis, based on files written by `experiment.py`.
# 
# `experiment.py` must be run first. 
#

import stats
import pandas.io.excel
import matplotlib.pyplot as plt

miceRAW = pandas.io.excel.read_excel('diet-nutrient-breakdown.xlsx','mice')
consumed = miceRAW['Dry weight food eaten g/mouse/cage/d']
(Cx, Cy) = stats.ecdf(consumed)

plt.clf()
	
p1, = plt.plot(Cx,Cy)
plt.xlabel('dry weight food eaten g/mouse/cage/day')
plt.ylabel('Proportion')
plt.gca().grid(True); # turn on grid lines.
plt.savefig('dryWeightEaten', dpi=600)