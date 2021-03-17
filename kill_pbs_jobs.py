# Specialial script for killing jobs on the University of Sydney high performance computing facility.
import subprocess

for i in list(range(2400868, 2401074)):
  cmd = 'qdel {:d}'.format(i)
  print(cmd)
  subprocess.call(['qdel', '{:d}'.format(i)])
