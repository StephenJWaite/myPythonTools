#Script is user to perform transforms between 0 position contour and the corresponding spiral contour
import numpy as np
import matplotlib.pyplot as plt
import imageTools as IT
from fieldwork.field.tools import fitting_tools

plt.switch_backend('TkAgg')
print '\nRunning script registerZeroPositionToSpiral...'
workingDir = './../imageToolsTestSpace/'

filename='6948_02'

#import the xyzPositions for corresponding spiral contour
spiralContour=np.loadtxt(workingDir+filename+'/CorrespondingSpiral/spiralControlPoints')

#import zero position contour from results
zPos,numContours,contours=IT.readContourFile(workingDir+filename+'/Results.txt')

#Ditch it all but keep the first contour
deformContour=contours[0]

xOpt,transformedContour=fitting_tools.fitTranslation(deformContour, spiralContour, xtol=1e-5, maxfev=0, sample=None, verbose=0, outputErrors=0)

#plot for sanity check
fig=plt.figure()
ax = fig.add_subplot(111)
ax.scatter(spiralContour[:,0],spiralContour[:,1]*-1,facecolors='g', edgecolors='g')
ax.scatter(deformContour[:,0],deformContour[:,1]*-1,facecolors='r', edgecolors='r')
#ax.scatter(transformedContour[:,0],transformedContour[:,1]*-1,facecolors='b', edgecolors='b')
ax.axis('equal')
for i in range(np.shape(spiralContour)[0]):
	print spiralContour[i,0]
	ax.plot([spiralContour[i,0],deformContour[i,0]],[spiralContour[i,1]*-1,deformContour[i,1]*-1])
	#ax.plot([spiralContour[i,0],deformContour[i,0]],[spiralContour[i,1],deformContour[i,1]]*-1)

plt.show()