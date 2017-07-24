#Script createControlPoints is used to palce control points on a png imgage. 
#Inputs: Png File directory
#        Contour directory
#        Number of control points
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


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

basePlot = plt.imshow(imageFile)
plt.show()


print 'Finished createControlPoints'

