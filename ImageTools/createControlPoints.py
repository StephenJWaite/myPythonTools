#Script createControlPoints is used to palce control points on a png imgage. 
#Inputs: Png File directory
#        Contour directory
#        Number of control points
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import imageTools as IT


print 'Running createControlPoints'

#Set up working directories
workingDir = '/home/stephen/Documents/GitRepos/myPythonTools/imageToolsTestSpace/'
imageDir = '/home/stephen/Documents/GitRepos/myPythonTools/imageToolsTestSpace/'
contourDir = '/home/stephen/Documents/GitRepos/myPythonTools/imageToolsTestSpace/'

#Read in input files
imageFile = mpimg.imread(imageDir + '6948_02.png')
print '    input Image size:', np.shape(imageFile)
contourData = np.loadtxt(contourDir + 'Contors.txt')
print '    input contour size:', np.shape(contourData)

#Sort the contour data
sortedContours=IT.sortContourData(contourData)

basePlot = plt.figure()
ax = basePlot.add_subplot(111)
ax.imshow(imageFile,cmap='gray') #Note Greys gives amazing contracts, may this for segmentation???

targetContour=sortedContours[0]

#TEST CODE - sub sample contour data for easy workability
print 'targetContour.shape:',np.shape(targetContour)
targetContour=targetContour[0::3]
print 'subsampled targetContour.shape:',np.shape(targetContour)

ax.plot(targetContour[:,0],targetContour[:,1],'-r')
ax.axis('equal')

#Set the number of master control Points
MCP=4
#Set the snap sphereRadius
sR=40
#Create a blank array for master points
MPoints=np.zeros((MCP,3))
#For the user specified number of master points, loop through and select.
for i in range(MCP):
	#Get the user input
	userPosition=plt.ginput(1)
	MPoints[i,:2]= np.asarray(userPosition)
	ax.scatter(MPoints[i,0],MPoints[i,1], facecolors='g', edgecolors='g')
	MPoints[i,:2],index = IT.snapSphere(targetContour,MPoints[i,:],sR)
	#insert the new points into the contour array
	#targetContour=np.concatenate((targetContour[:index+1,:],np.asarray([xPt,yPt]).reshape(1,2),targetContour[index+1:,:]))
	targetContour=np.concatenate((targetContour[:index+1,:],MPoints[i,:2].reshape(1,2),targetContour[index+1:,:]))
	MPoints[i,2]=index+1
	ax.scatter(MPoints[i,0],MPoints[i,1], facecolors='y', edgecolors='y')

#Seed the slave control points between the masters with equidistance spacing.
print 'Master Points\n',MPoints
distVector=IT.calculateDistanceVector(targetContour,MPoints)
print 'Dist Vector\n',distVector
seedNumber=2
controlPoints=IT.seedSlavePoints(distVector,contourData,MPoints,seedNumber)
plt.show()


print 'Finished createControlPoints'

