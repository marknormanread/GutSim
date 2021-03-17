'''
Created on 10/10/2014

@author: markread
'''
from random import shuffle
import random

def shuffle_subset(l, start, end):
	'''
	Start and end are the indexes between which to shuffle. Start is inclusive, end is exclusive (1 greater than the
	index comprising the end of the slice.
	'''
	copy = l[start:end]
	print copy
	shuffle(copy)
#	l[start:end] = copy


a = range(205)

intervals = 10
shuffleSize = int(len(a) / (intervals + 1))
print 'shuffleSize = ' + str(shuffleSize)
startPeriodShuffle = random.randrange(shuffleSize)
shuffle_subset(a, 0, shuffleSize)
shuffle_subset(a, -shuffleSize, len(a))
i = 0
while True:
	start = (i * shuffleSize) + startPeriodShuffle
	end = start + shuffleSize
	print 'start ' + str(start)
	print 'end  ' + str(end)
	if end > len(a):
		break
	shuffle_subset(a, start, end)
	i += 1

print a