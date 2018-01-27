from pyevtk.hl import pointsToVTK
import numpy as np

workingDir='/home/stephen/OpenFOAM/Simulations2/MotionCheck/'

#read in some points
points=np.loadtxt(workingDir+'targetPoints_5.txt',delimiter=',')

#write vtk
#pointsToVTK("./test_vtk",points[:,0],points[:,1],points[:,2],data={"cat" : points[:,2]})

x = np.arange(1.0,10.0,0.1)
y = np.arange(1.0,10.0,0.1)
z = np.arange(1.0,10.0,0.1)
print z
print np.shape(x),np.shape(y),np.shape(z)

print np.shape(points[:,0]),np.shape(points[:,1]),np.shape(points[:,2])
print points[1:5,:]

pointsToVTK("./line_points", x, y, z, data = {"elev" : z.T})

print x.flags
print '--------------------'
print points.flags
points=np.asfortranarray(points)
print points.flags
pointsToVTK(workingDir+"test_vtk",points[1:5,0],points[1:5,1],points[1:5,1],data={"cat" : points[1:5,1]})