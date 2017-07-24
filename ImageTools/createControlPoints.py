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
workingDir = '/home/stephen/Documents/Python/imageToolsTestSpace/'
imageDir = '/home/stephen/Documents/Python/imageToolsTestSpace/'
contourDir = '/home/stephen/Documents/Python/imageToolsTestSpace/'

#Read in input files
imageFile = mpimg.imread(imageDir + '6948_02.png')
print '    input Image size:', np.shape(imageFile)
contourData = np.loadtxt(contourDir + 'Contors.txt')
print '    input contour size:', np.shape(contourData)

basePlot = plt.figure()
fuck = basePlot.add_subplot(111)
fuck.imshow(imageFile)

#rezero the slide numbers imageFile[:,2]
startSlide = np.min(contourData[:,2])
print 'Starting slide:',startSlide

#loop through and plot the first contour
for i in range(len(contourData)):
	if contourData[i,2] <= startSlide:
		fuck.scatter(contourData[i,0]+5,contourData[i,1]+5)

fuck.axis('equal')

points = IT.pickPoints(3,basePlot)

plt.show()


print 'Finished createControlPoints'

