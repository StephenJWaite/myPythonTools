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

workingDir = '/home/stephen/OpenFOAM/PhDSimulations/Motion2/ContractionContours/'

#we will have a control file that has all of the contour groups we are planning to use
activeContours=open(workingDir+'BiContourList').readlines()

#read in the contraction map
CMap=np.loadtxt(workingDir+'contractionMap')

#set up infromation arrays
contourData=[]
contourInformation=[]
#loop through activeContours
for line in activeContours:
	print 'Opening contour:',line.rstrip()
	zPos,numContours,contours=IT.readContourFile(workingDir+line.rstrip()+'/Results.txt')
	#store the zPos and numContous in contourInformation
	contourInformation=contourInformation+[[zPos,numContours]]
	contourData=contourData+[contours]

contractionMatrix,endtime=IT.createContractionMatrix(activeContours,contourInformation,contourData,CMap)

print np.shape(contractionMatrix)
print np.shape(contourInformation)
print contourInformation[0][0]
#print contourInformation[1][0]
#print contourInformation[2][0]
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
	timePlot.set_xlim(-140,75)
	timePlot.set_ylim(-140,150)
	timePlot.view_init(elev=-24, azim=-163)
	for j in range(np.shape(contractionMatrix)[1]):
		timePlot.plot(contractionMatrix[i][j][:,0],contractionMatrix[i][j][:,1],np.float(contourInformation[j][0]),'-or')

	plt.pause(0.2)

fig3=plt.figure(3)
contourPlot=fig3.add_subplot(111)
for i in range(np.shape(contractionMatrix)[1]):
	for j in range(endtime):
		contourPlot.clear()
		contourPlot.set_xlim(-140,75)
		contourPlot.set_ylim(-140,150)
		contourPlot.plot(contractionMatrix[j][i][:,0],contractionMatrix[j][i][:,1]*-1,np.float(contourInformation[i][0]),'-or')
		plt.pause(0.2)

	


plt.ioff()

#Lets now look in 3D
#fig3=plt.figure(3)
#ax3D=fig3.add_subplot(111, projection='3d')
#for imLoop in range(3,endtime-3,2):
#    ax3D.plot(contractionMatrix[imLoop][0][:,0],contractionMatrix[imLoop][0][:,1],imLoop,'-or')

#plt.show()
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
