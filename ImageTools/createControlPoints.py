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

basePlot = plt.figure()
ax = basePlot.add_subplot(111)
ax.imshow(imageFile,cmap='gray') #Note Greys gives amazing contracts, may this for segmentation???

#rezero the slide numbers imageFile[:,2]
startSlide = np.min(contourData[:,2])
print 'Starting slide:',startSlide

#loop through and plot the first contour
targetContour=contourData[0,0:2].reshape(1,2)
for i in range(len(contourData)):
	if contourData[i,2] <= startSlide:
		targetContour=np.append(targetContour,contourData[i,0:2].reshape(1,2),axis=0)

print 'targetContour.shape:', np.shape(targetContour)

#TEST CODE - sub sample contour data for easy workability
targetContour=targetContour[0::4]
print 'subsampled targetContour.shape: ',np.shape(targetContour)

ax.plot(targetContour[:,0],targetContour[:,1],'-or')
ax.axis('equal')

#Set the number of master control Points
MCP=4
#Set the snap sphereRadius
sR=40
#Create a blank array for master points
MPoints=np.zeros((MCP,2))
#set up the new array 
for i in range(MCP):
	userPosition=plt.ginput(1)
	print 'Point group ',i ,':', userPosition[0]
	MPoints[i,:]= np.asarray(userPosition)
	ax.scatter(MPoints[i,0],MPoints[i,1], facecolors='g', edgecolors='g')
	xPt,yPt, index = IT.snapSphere(targetContour,MPoints[i,:],sR)
	#insert the new points into the contour array
	temp=np.zeros((len(targetContour),2))
	np.copyto(temp,targetContour)
	print '+++++++++++++++++++++++++++++++++++++++'
	print 'array size: ',np.shape(targetContour)
	targetContour=np.concatenate((temp[:index+1,:],np.asarray([xPt,yPt]).reshape(1,2),temp[index+1:,:]))
	print 'New array size: ',np.shape(targetContour)
	print 'inserted index: ', index
	print 'array values:\n',targetContour
	print '+++++++++++++++++++++++++++++++++++++++'
	ax.scatter(xPt,yPt, facecolors='y', edgecolors='y')

#Sanity plot the new contour to look for twist

plt.show()


print 'Finished createControlPoints'

