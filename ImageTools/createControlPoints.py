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
workingDir = '/home/stephen/Documents/Python/myPythonTools/imageToolsTestSpace/'
imageDir = '/home/stephen/Documents/Python/myPythonTools/imageToolsTestSpace/'
contourDir = '/home/stephen/Documents/Python/myPythonTools/imageToolsTestSpace/'

#Read in input files
imageFile = mpimg.imread(imageDir + '6948_02.png')
print '    input Image size:', np.shape(imageFile)
contourData = np.loadtxt(contourDir + 'Contors.txt')
print '    input contour size:', np.shape(contourData)

basePlot = plt.figure()
ax = basePlot.add_subplot(111)
ax.imshow(imageFile)

#rezero the slide numbers imageFile[:,2]
startSlide = np.min(contourData[:,2])
print 'Starting slide:',startSlide

#loop through and plot the first contour
firstContour=contourData[0,0:2].reshape(1,2)
for i in range(len(contourData)):
	if contourData[i,2] <= startSlide:
		firstContour=np.append(firstContour,contourData[i,0:2].reshape(1,2),axis=0)

print 'firstContour.shape:', np.shape(firstContour)
ax.plot(firstContour[:,0]+5,firstContour[:,1]+5)

ax.axis('equal')

#Set the number of master control Points
MCP=4
#Set the snap sphereRadius
sR=20
MPoints=np.zeros((MCP,2))
for i in range(MCP):
	cat=plt.ginput(1)
	print 'Point group ',i ,':', cat[0]
	MPoints[i,:]= np.asarray(cat)
	ax.scatter(MPoints[i,0],MPoints[i,1], facecolors='none', edgecolors='g')
	IT.snapSphere(firstContour,MPoints[i,:],sR)

plt.show()


print 'Finished createControlPoints'

