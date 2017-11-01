#Script createControlPoints is used to palce control points on a png imgage.
#Inputs: Png File directory
#        Contour directory
#        Number of control points
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.image as mpimg
import imageTools as IT
import time

plt.switch_backend('TkAgg')

print 'Running createControlPoints'

#Set up working directories
location='BlindSacs'
fileName='6905a_02'
contraction='Contraction2'
region='Dorsal'

workingDir = '/home/stephen/Documents/PhD/Data/Motility/SortedByLocation/'+location+'/'+fileName+'/'
imageDir = '/home/stephen/Documents/PhD/Data/Motility/SortedByLocation/'+location+'/'+fileName+'/Slices/'+contraction+'/'
contourDir = '/home/stephen/Documents/PhD/Data/Motility/SortedByLocation/'+location+'/'+fileName+'/'



#Read in input files
contourData = np.loadtxt(contourDir + contraction+region)
print '    input contour size:', np.shape(contourData)
#Sort the contour data
sortedContours=IT.sortContourData(contourData)
ControlPointsArray=[0]*np.shape(sortedContours)[0]

basePlot = plt.figure()

 #Note Greys gives amazing contraction contrast, mayne this for segmentation???

#Set the number of master control Points
MCP=1 #4
#Set the snap sphereRadius
sR=40
#set the number of slave control points per master point
seedNumber=35#8 

#Begin loop
for imLoop in range(np.shape(sortedContours)[0]):
    targetContour=sortedContours[imLoop]

    ax = basePlot.add_subplot(111)
    imageFile = mpimg.imread(imageDir + str(imLoop) +'.png')
    print '    input Image size:', np.shape(imageFile)
    ax.imshow(imageFile,cmap='gray')
    plt.draw()
	#TEST CODE - sub sample contour data for easy workability
    #print '\ttargetContour.shape:',np.shape(targetContour)
    #targetContour=targetContour[0::5]
    #print '\tsubsampled targetContour.shape:',np.shape(targetContour)

    print np.size(targetContour), targetContour
    if np.size(targetContour)==1:
        print 'BLARP'
        ControlPointsArray[imLoop]=targetContour

    else:
        ax.plot(targetContour[:,0],targetContour[:,1],'-xr')
        ax.axis('equal')

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
            #dirty shift fix, should recode this whole section to first determine
            #where the M points should all go, and then place them in order after
            #wards, but since num MP is low, its not gonna slow stuff down.
            for j in range(i):
                if MPoints[i,2]<MPoints[j,2]:
                    MPoints[j,2]=MPoints[j,2]+1

            ax.scatter(MPoints[i,0],MPoints[i,1], facecolors='y', edgecolors='y')

        #Seed the slave control points between the masters with equidistance spacing.
        print 'Master Points\n',MPoints
        distVector=IT.calculateDistanceVector(targetContour,MPoints)
        print 'Dist Vector\n',distVector

        controlPoints=IT.seedSlavePoints(distVector,targetContour,MPoints,seedNumber,ax)
        ax.plot(controlPoints[:,0],controlPoints[:,1],'ob')
        ControlPointsArray[imLoop]=controlPoints
    
    plt.draw()
    time.sleep(0.5)
    basePlot.clf()
    

#Lets now look in 3D
#fig2=plt.figure(2)
#ax3D=fig2.add_subplot(111, projection='3d')
#for imLoop in range(np.shape(sortedContours)[0]):
#    ax3D.plot(ControlPointsArray[imLoop][:,0],ControlPointsArray[imLoop][:,1],imLoop,'-or')

plt.show()

#Write outfile
#Set the z position
#Set the total time (number of contours)
with file(workingDir+fileName+'_'+contraction+region,'w') as outfile:
    outfile.write('#Z position: {0}\n'.format(5))
    outfile.write('#Contraction length: {0}\n'.format(np.shape(sortedContours)[0]))
    for i in range(np.shape(sortedContours)[0]):

        if np.size(ControlPointsArray[i])==1:
            print 'RARP'
            cX=cY=0
            print 'Centroid pos', cX,cY
            CS=0
            print 'Centroid Size', CS
            Area=0
            print 'Area:', Area
            outfile.write('#Contour: {0} Centroid: [{1},{2}] Centroid Size: {3} Area: {4}\n'.format(i,cX,cY,CS,Area))
        else:
            cX,cY=IT.calculateCentroidPosition(ControlPointsArray[i][:,0],ControlPointsArray[i][:,1])
            print 'Centroid pos', cX,cY
            CS=IT.calculateCentroidSize(ControlPointsArray[i][:,0],ControlPointsArray[i][:,1])
            print 'Centroid Size', CS
            Area=IT.calculateClosedContourArea(ControlPointsArray[i][:,0],ControlPointsArray[i][:,1])
            print 'Area:', Area
            outfile.write('#Contour: {0} Centroid: [{1},{2}] Centroid Size: {3} Area: {4}\n'.format(i,cX,cY,CS,Area))
            np.savetxt(outfile,ControlPointsArray[i])

#with file(workingDir+'Results.txt','w') as outfile:
	#outfile.write('#Z position: {0}\n'.format(5))
	#outfile.write('#Contraction length: {0}\n'.format(np.shape(sortedContours)[0]))   
    #for i in range(np.shape(sortedContours)[0]):
        #print i
        #cX,cY=IT.calculateCentroidPosition(ControlPointsArray[i][:,0],ControlPointsArray[i][:,1])
        #print cX,cY
        #CentroidSize=IT.calculateCentroidSize(ControlPointsArray[i][:,0],ControlPointsArray[i][:,1])
        #print CentroidSize
        #Area=IT.calculateClosedContourArea(ControlPointsArray[i][:,0],ControlPointsArray[i][:,1])
        #print Area
        #outfile.write('#Contour: {0} Centroid: [{0},{0}] Centroid Size: {0} Area: {0}'.format(i,cX,cY,CentroidSize,Area))
        #np.savetxt(outfile,ControlPointsArray[i])

print 'Finished createControlPoints'