
'''
Created on 09/09/2014

@author: markread



'''
import sys
import numpy as np
import scipy.stats as ss
import os.path

import compileMedians
import stats

### used for writing response analysis data as an xml file ###
from xml.etree.ElementTree import Element, SubElement
from xml.etree import ElementTree
from xml.dom import minidom



def setupResponse(name, parent, ksd, ksp, a):
	'''
	Used in constructing XML trees.
	Creates tags representing response data. Includes data for KS tests (p and d) and A test.
	'''
	child = SubElement(parent, name)
	childKSp = SubElement(child,'ks-p')
	childKSp.text = str(ksp)
	childKSd = SubElement(child,'ks-d')
	childKSd.text = str(ksd)
	childA = SubElement(child,'A')
	childA.text = str(a)


def prettify(elem):
	'''
	Return a pretty-printed XML string for the Element.

	Taken from: http://pymotw.com/2/xml/etree/ElementTree/create.html on 19/08/2014
	'''
	rough_string = ElementTree.tostring(elem, 'utf-8')
	reparsed = minidom.parseString(rough_string)
	return reparsed.toprettyxml(indent="  ")




############################ PROGRAM STARTS HERE ############################

directory = '.'
directory = '../results/aleatory/5'
if len(sys.argv) > 1:
	directory = sys.argv[1]


paramVals = range(20)

aTestBoundsMin = 0.29;
aTestBoundsMax = 0.71;


# read in files containing data. Each single file can contain several lines describing replicate experiments.
data = None
headerString = None		# the header at top of results file.
for i,p in enumerate(paramVals):
	inFile = directory + '/' + str(p) + '/all-mouse'

	if not os.path.isfile(inFile):
		# medians data not compiled, compile it now.
		compileMedians.compile(str(p))

	if headerString is None:
		with open(inFile,'r') as f1:
			headerString = f1.readline()

	d = np.genfromtxt(inFile, dtype=float, delimiter=',', skip_header=1)
	if data is None:
		dim = [len(paramVals)]
		dim.extend(d.shape)
		data = np.zeros(dim)
	data[i,:,:] = d

# split into array of column names. Drop the '#' at start, and remove '\n' at end.
header = headerString[1:-1].split(',')

analXML = Element('analysis')

default = 0		# the index of the defailt parameter set in the data array.
As = []			# store A test scores in here. Order same as paramVals.
for i,p in enumerate(paramVals):
	paramXML = SubElement(analXML,'parameter')
	paramValXML = SubElement(paramXML,'value')
	paramValXML.text = str(p)

	for j,h in enumerate(header):
		safeH = h.replace('/','-')			# replace illegal characters in XML tag to something else.
		defaultDistro = data[default,:,j]
		paramDistro = data[i,:,j]
		a = stats.Atest(defaultDistro,paramDistro)
		(ksd, ksp) = ss.ks_2samp(defaultDistro,paramDistro)
		setupResponse(name=safeH,parent=paramXML,ksd=ksd,ksp=ksp,a=a)	# write response data to XML file.

fout = open(directory + '/' + 'robustnessAnalysisResults.xml','w')
fout.write(prettify(analXML))
fout.close()
print prettify(analXML)


# plot graph of responses
responseXPaths = [	'parameter/mouse-_ID/A',	# this works.
					'./parameter//relBcf/A',	# paths to identify response values in the XML file.
					'.//relBcm/A',
					'.//relBdf/A',
					'.//relBdf/A',
					'.//relBmm/A',
				]

for path in responseXPaths:
	print analXML.findall(path)
