import cProfile
import re
import gutsim.experiment
import pstats

cProfile.run('gutsim.experiment.profile()', 'profile_stats_raw')
stream = open('profile_stats','w')
p = pstats.Stats('profile_stats_raw', stream=stream)
p.strip_dirs().sort_stats('tottime').print_stats()


