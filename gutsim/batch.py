# '''
# Created on 01/05/2014

# @author: markread


# This module represents a batch culture of bacteria
# '''
# import gutsim.vessel
# import gutsim.bacteria
# from gutsim.nutrients import Nutrients
# from random import shuffle
# import math


# class Batch(gutsim.vessel.Vessel):
# 	def __init__(self, initBcf, initBcm, initBdf, initBdm, initBmm, Cc, Cd, Cm, Nf, Nm, *args, **kwargs):
# 		super(Batch, self).__init__(*args, **kwargs)
# 		self._initBcf = initBcf
# 		self._initBcm = initBcm
# 		self._initBdf = initBdf
# 		self._initBdm = initBdm
# 		self._initBmm = initBmm
# 		self._Cc = Cc
# 		self._Cd = Cd
# 		self._Cm = Cm
# 		self._Nf = Nf
# 		self._Nm = Nm
# 		self._nutrients = Nutrients(Cw=Cc, Cd=Cd, Cm=Cm, Nf=Nf, Nm=Nm)

# 	def _uniform_death(self, cumDied, timestep):
# 		'''
# 		Labelled uniform beacuse all bacteria have an equal chance of dying under this regime. There is no notions that
# 		some bacteria may be more resistant to death than others.
# 		'''
# 		# Logarithmic decay. dN/dt = -tao . N
# 		# transit time represents mean lifetime (ML).
# 		# hence, death rate (tao) = 1 / (ML).
# 		meanLifetime = 18.0
# 		deathRate = 1.0 / meanLifetime
# 		# keep a cumulative count, because after scaling by timeslice the quantity to be removed may be < 1.
# 		cumDied += len(self._bacteria) * deathRate * timestep	# scale by timestep
# 		remove = math.floor(cumDied)		# can only remove whole numbers of bacteria
# 		cumDied -= remove
# 		for _ in range(int(remove)):
# 			del(self._bacteria[-1])	# delete off end of list
# 		return cumDied

# 	def execute(self, endtime, timestep=0.1):
# 		self.logger = gutsim.vessel.Logger(vessel=self, prefix='batch/', logPeriod=10,bacterGraphs=False,graphPath='batch/bacteria/')
# 		# initialise bacteria in the simulation
# 		for _ in range(self._initBcf):
# 			self._bacteria.append(gutsim.bacteria.Guild_rf_f(host=self, timestep=timestep, carbon=10.0, nitrogen=10.0))
# 		for _ in range(self._initBcm):
# 			self._bacteria.append(gutsim.bacteria.Bcm(host=self, timestep=timestep, carbon=10.0, nitrogen=10.0))
# 		for _ in range(self._initBdf):
# 			self._bacteria.append(gutsim.bacteria.Bdf(host=self, timestep=timestep, carbon=10.0, nitrogen=10.0))
# 		for _ in range(self._initBdm):
# 			self._bacteria.append(gutsim.bacteria.Bdm(host=self, timestep=timestep, carbon=10.0, nitrogen=10.0))
# 		for _ in range(self._initBmm):
# 			self._bacteria.append(gutsim.bacteria.Bmm(host=self, timestep=timestep, carbon=10.0, nitrogen=10.0))
# 		currentTime = 0.0
# 		self.logger.log(self._bacteria, self._nutrients, currentTime)
# 		while currentTime < endtime:
# 			currentTime += timestep
# 			if True: print('ct = ' + str(currentTime) + ' ' + self.logger.string_state())

# 			if str(currentTime) == '30000.0':
# 				print("now")
# 				self._nutrients.add_quantities(Cc=self._Cc, Cd=self._Cd, Cm=self._Cm, Nf=self._Nf, Nm=self._Nm)

# 			shuffle(self._bacteria)
# 			for b in self._bacteria:
# 				consumed = b.step(self._nutrients)
# 				self._nutrients.subtract(consumed)

# 			self._bacteria.extend(self._newBacteria)	# register the new bacteria.
# 			self._newBacteria = []						# clear the new bacteria register for next step.
# 			# use set difference to delete expired bacteria.
# 			#self._bacteria = list(set(self._bacteria) - set(self._deadBacteria))
# 			self._deadBacteria = []						# reset for next step

# 			# perform logging of bacterial numbers.
# 			self.logger.log(self._bacteria, self._nutrients, currentTime)



# 	def stringBacteriaState(self):
# 		for b in self._bacteria:
# 			print(str(b))


# def main():
# 	Cc = 0.1
# 	Cd = 0.1
# 	Cm = 0.1
# 	Nf = 0.1
# 	Nm = 0.1

# 	batch = Batch(initBcf=100, initBcm=0, initBdf=0, initBdm=0, initBmm=0,
# 				  Cc=Cc, Cd=Cd, Cm=Cm, Nf=Nf, Nm=Nm)
# 	batch.execute(endtime=100, timestep=0.1)			# 7 days = 168 hours
# 	#batch.logger.writeMouseToFile('batch/abundances')
# 	batch.logger.graph_global_dynamics('batch/timeseries.png')
# 	batch.logger.graph_guild_states(prefix='batch/')

# if  __name__ =='__main__':main()