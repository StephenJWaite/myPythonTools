#Script createControlPoints is used to palce control points on a png imgage.
#Inputs: Png File directory
#        Contour directory
#        Number of control points
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pltLab
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.image as mpimg
import imageTools as IT

plt.switch_backend('TkAgg')

print 'Running createControlPoints for single case'

#Set up working directories
workingDir = './../imageToolsTestSpace/6948_02/CorrespondingSpiral/'
imageDir = './../imageToolsTestSpace/6948_02/CorrespondingSpiral/'
contourDir = './../imageToolsTestSpace/6948_02/CorrespondingSpiral/'

#Read in input files
contourData = np.loadtxt(contourDir + 'Contour.txt')
print '    input contour size:', np.shape(contourData)
#Sort the contour data

basePlot = plt.figure()

#Set the number of master control Points
MCP=4
#Set the snap sphereRadius
sR=40
#set the number of slave control points per master point
seedNumber=8
  
ax = basePlot.add_subplot(111)
imageFile = mpimg.imread(imageDir + '524.png')
print '    input Image size:', np.shape(imageFile)
ax.imshow(imageFile,cmap='gray')

ax.plot(contourData[:,0],contourData[:,1],'-xr')
ax.axis('equal')

#Create a blank array for master points
MPoints=np.zeros((MCP,3))

#For the user specified number of master points, loop through and select.
for i in range(MCP):
    #Get the user input
    userPosition=plt.ginput(1)
    MPoints[i,:2]= np.asarray(userPosition)
    ax.scatter(MPoints[i,0],MPoints[i,1], facecolors='g', edgecolors='g')
    MPoints[i,:2],index = IT.snapSphere(contourData,MPoints[i,:],sR)
    #insert the new points into the contour array
    #contourData=np.concatenate((contourData[:index+1,:],np.asarray([xPt,yPt]).reshape(1,2),contourData[index+1:,:]))
    contourData=np.concatenate((contourData[:index+1,:],MPoints[i,:2].reshape(1,2),contourData[index+1:,:]))
    MPoints[i,2]=index+1
    #dirty shift fix, should recode this whole section to first determine
    #where the M points should all go, and then place them in order after
    #wards, but since num MP is low, its not gonna slow stuff down.
    for j in range(i):
        if MPoints[i,2]<MPoints[j,2]:
            MPoints[j,2]=MPoints[j,2]+1

    ax.scatter(MPoints[i,0],MPoints[i,1], facecolors='y', edgecolors='y')

#Seed the slave control points between the masters with equidistance spacing.
print 'Master Points\n',MPoints
distVector=IT.calculateDistanceVector(contourData,MPoints)
print 'Dist Vector\n',distVector

controlPoints=IT.seedSlavePoints(distVector,contourData,MPoints,seedNumber,ax)
ax.plot(controlPoints[:,0],controlPoints[:,1],'ob')
plt.draw()



#Calculate the centroid size
cx,cy=IT.calculateCentroidPosition(controlPoints[:,0],controlPoints[:,1])
ax.plot(cx,cy,'or')
CS=IT.calculateCentroidSize(controlPoints[:,0],controlPoints[:,1])
print 'Centroid Size is',CS
#Calculate the area of the contour
Area1=IT.calculateClosedContourArea(controlPoints[:,0],controlPoints[:,1])
print Area1/1e6
Area2=IT.PolygonArea(np.hstack([controlPoints[:,0],controlPoints[0,0]]),np.hstack([controlPoints[:,1],controlPoints[0,1]]))
print Area2/1e6

#Write outfile
np.savetxt(workingDir+'spiralControlPoints',controlPoints)


print 'Finished createControlPoints'

plt.show()