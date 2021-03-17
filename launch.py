"""
Created on 16/09/2014

Entry point into the simulation.

@author: Mark N. Read
"""
import sys
if '-longitudinal' in sys.argv:
    import gutsim.experiment_longitudinal
    gutsim.experiment_longitudinal.main()
else:
    import gutsim.experiment
    gutsim.experiment.main()