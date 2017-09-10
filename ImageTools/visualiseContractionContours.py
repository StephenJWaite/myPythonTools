#Script visualiseContractionContours.py is used as a sanity check for contour information to see how it will be applied
#within Anthony.
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation
import imageTools as IT
import time

plt.switch_backend('TkAgg')

print '\nRunning script visualiseContractionContours...'

workingDir = './../imageToolsTestSpace/'

#we will have a control file that has all of the contour groups we are planning to use
activeContours=open(workingDir+'ContourList').readlines()
contourData=[]
contourInformation=[]

#loop through activeContours
for line in activeContours:
	print 'Opening contour:',line.rstrip()
	zPos,numContours,contours=IT.readContourFile(workingDir+line.rstrip()+'/Results.txt')
	#store the zPos and numContous in contourInformation
	contourInformation=contourInformation+[[zPos,numContours]]
	contourData=contourData+[contours]

#Print out some sanity statistics
print 'Contour information:'
for i in range(np.size(activeContours)):
	print '\t',activeContours[i].rstrip(),'zPos:',contourInformation[i][0],'number Contours',contourInformation[i][1]

#Loop through and plot in 3D each contours slices at the correct zPos


#So the current plan is, to bit a time range (roughly 2 minutes).
#we create a nxd array where n is the number of slice planes we are using, and d is seconds.
#using our contraciton map from literature, we will place active contours in the correct
#times, and duplicate their first and last slices to rill out time on either side.

#read in the contraction map
CMap=np.loadtxt(workingDir+'contractionMap')

#will need to calculate a suitable end time, as some contours contract for a while.
#This will be dependent on when they start in the contraciton map. 
#For now we will take the largest contraciton time and tack it onto the end time in the 
#CMAP
endtime=np.int(CMap[-1,1]+max(np.asarray(contourInformation)[:,1].astype(float)))+1
print 'End time is:',endtime


#set up the contractionMatrix
#contractionMatrix=[[0]*np.shape(activeContours)[0]]*endtime <-never do this
contractionMatrix=[[0]*np.shape(activeContours)[0] for i in range(endtime)]
#
print 'contractionMatrix shape:',np.shape(contractionMatrix)
MapFunc=interp1d(CMap[:,0],CMap[:,1],kind='cubic')

for i in range(np.shape(activeContours)[0]):
	spoint=np.int(MapFunc(np.float(contourInformation[i][0])))
	print 'Processing slice:',i,'at',contourInformation[i][0]
	print '\tStart Points:',spoint, 'End Point:',spoint+int(contourInformation[i][1])
	print 'looping from',0,'to',spoint
	for j in range(spoint):
		#print 'i:',i,'j:',j,'slice:',0
		contractionMatrix[j][i]=contourData[i][0]

	print 'looping from',spoint,'to',np.int(spoint+int(contourInformation[i][1]))
	for j in range(spoint,np.int(spoint+int(contourInformation[i][1]))):
		#print 'i:',i,'j:',j,'slice',(j-spoint)
		contractionMatrix[j][i]=contourData[i][j-spoint]

	print 'looping from',np.int(spoint+int(contourInformation[i][1])),'to',endtime
	for j in range(np.int(spoint+int(contourInformation[i][1])),endtime):
		#print 'i:',i,'j:',j,'slice: end'
		contractionMatrix[j][i]=contourData[i][-1]

#fig=plt.figure()
#ax3D=fig.add_subplot(111, projection='3d')

#for i in range(endtime):
#	ax3D.plot(contractionMatrix[i][0][:,0],contractionMatrix[i][0][:,1],i,'-or')

#Some methods for temporal ploting#

#####################################################################
#Brute force with plot options, this totaly works though #Baller
fig2=plt.figure(2)
timePlot=fig2.add_subplot(111,projection='3d')
timePlot.set_xlim(50,350)
timePlot.set_ylim(100,400)
plt.ion()

for i in range(endtime):
	timePlot.clear()
	timePlot.set_xlim(50,350)
	timePlot.set_ylim(100,400)
	timePlot.view_init(elev=-81+i, azim=-113+i)
	timePlot.plot(contractionMatrix[i][0][:,0],contractionMatrix[i][0][:,1],5,'-or')
	timePlot.plot(contractionMatrix[i][1][:,0],contractionMatrix[i][1][:,1],1,'-or')
	plt.pause(0.05)

plt.ioff()
plt.show()
#######################################################################

######################################################################
#Using animation functionality in matplotlib
#testFig=plt.figure(2)
#ax = plt.axes(xlim=(0, 2), ylim=(-2, 2))
#line,=ax.plot([],[],'r-')

#def update(i,):
#	line.set_data(contractionMatrix[i][0][:,0],contractionMatrix[i][0][:,1])
#	return line

#anim=matplotlib.animation.FuncAnimation(testFig,update,frames=10,repeat=True)
#plt.show()
########################################################################



















