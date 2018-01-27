#this is getting crazy, so this script is to deform the geometry to some % of contracted, and mesh the geometry
import os,sys
sys.path.insert(0,'./../../lib')
import numpy as np
import scipy
import imageTools as IT
from pyevtk.hl import pointsToVTK
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#HMF depedancies
from scipy.spatial import cKDTree
from fieldwork.field.tools import fitting_tools
from fieldwork.field import geometric_field
from fieldwork.field import geometric_field_fitter as GFF
from gias.common import fieldvi, transform3D
from gias.common import alignment_fitting as af
#Set up the working dirs
workingDir = '/home/stephen/OpenFOAM/Simulations2/MotionCheck/RemeshAttempt/Grouped/'
workingDirOF = '/home/stephen/OpenFOAM/Simulations2/MotionCheck/RemeshAttempt/RumenRoofInlet'

#we will have a control file that has all of the contour groups we are planning to use
activeContours=open(workingDir+'ContourList').readlines()

#read in the contracted position
contractionPositions=open(workingDir+'MeshPositions').readlines()
source_points=np.empty((0,3))
target_points=np.empty((0,3))

for line in contractionPositions:
	#read in the apropraite file
	print line.rstrip().split()[0]
	#temp=np.loadtxt(workingDir+line.rstrip().split()[0]+'/Results.txt')
	zPos,numContours,contours=IT.readContourFile(workingDir+line.rstrip().split()[0]+'/Results.txt')
	temp=np.ones((np.shape(contours)[1],1))*np.float(zPos)
	#print np.append(contours[int(line.rstrip().split()[1])],temp,axis=1)
	source_points=np.vstack((source_points,np.append(contours[0],temp,axis=1)))
	target_points=np.vstack((target_points,np.append(contours[int(line.rstrip().split()[1])],temp,axis=1)))

print source_points.shape

#Add in some fixed contours
fixedContour1=np.loadtxt(workingDir+'FixedContours/C1')
source_points=np.vstack((source_points,fixedContour1))
target_points=np.vstack((target_points,fixedContour1))

scalar=np.ones(np.shape(target_points)[0])
#print 'SCALAR SHAPE!!!!!',scalar
points=np.asfortranarray(target_points)
pointsToVTK(workingDirOF+"/targetPoints",points[:,0],points[:,1],points[:,2],data={"cat" : scalar})


#plot results
fig1=plt.figure(2)
Plot=fig1.add_subplot(111,projection='3d')
Plot.plot(source_points[:,0],source_points[:,1],source_points[:,2],'or')
Plot.plot(target_points[:,0],target_points[:,1],target_points[:,2],'ob')

#now produce a stl for meshing
# fititng parameters for host mesh fitting
host_mesh_pad = 20.0 # host mesh padding around slave points
host_elem_type = 'quad555' # quadrilateral cubic host elements
host_elems = [1,1,1] # a single element host mesh
maxit = 20
sobd = [4,4,4]
sobw = 1e-10
xtol = 1e-12

#import the ply points
wallPoints = np.loadtxt(workingDirOF+'/Wall.xyz')
inletPoints=np.loadtxt(workingDirOF+'/Inlet.xyz')

# make host mesh
host_mesh = GFF.makeHostMeshMulti(
            	wallPoints.T,
                host_mesh_pad,
                host_elem_type,
                host_elems,
                )

#embedding each patch For passive deform
wallPoints_xi = host_mesh.find_closest_material_points(
                            wallPoints,
                            initGD=[50,50,50],
                            verbose=True,
                            )[0]

wallPoints_passive = geometric_field.makeGeometricFieldEvaluatorSparse(
								host_mesh, [1,1],
    							matPoints=wallPoints_xi,
                                )

inletPoints_xi = host_mesh.find_closest_material_points(
                            inletPoints,
                            initGD=[50,50,50],
                            verbose=True,
                            )[0]

inletPoints_passive = geometric_field.makeGeometricFieldEvaluatorSparse(
								host_mesh, [1,1],
    							matPoints=inletPoints_xi,
                                )

print 'calculating slave xi...'
slave_xi = host_mesh.find_closest_material_points(
					source_points,
					initGD=[100,100,100],
					verbose=True,
					)[0]

target_tree = cKDTree(target_points)

def my_slave_func(v2):
		v1=target_points
		d=np.sqrt((v1[:,0]-v2[:,0])**2 + (v1[:,1]-v2[:,1])**2 + (v1[:,2]-v2[:,2])**2)
		return d

slave_func = my_slave_func

# host mesh fit
host_x_opt, source_points_fitting_hmf,\
slave_xi, rmse_hmf = fitting_tools.hostMeshFitPoints(
                host_mesh,
                source_points,
                slave_func,
                slave_xi=slave_xi,
                max_it=maxit,
                sob_d=sobd,
                sob_w=sobw,
                verbose=True,
                xtol=xtol
                )
	
# evaluate the new positions of the passive source points
wallPoints_hmf = wallPoints_passive(host_x_opt).T
inletPoints_hmf = inletPoints_passive(host_x_opt).T

scalar=np.ones(np.shape(target_points)[0])
#print 'SCALAR SHAPE!!!!!',scalar
points=np.asfortranarray(target_points)
pointsToVTK(workingDirOF+"/targetPoints",points[:,0],points[:,1],points[:,2],data={"cat" : scalar})
points=np.asfortranarray(source_points_fitting_hmf)
pointsToVTK(workingDirOF+"/sourcePoints",points[:,0],points[:,1],points[:,2],data={"cat" : scalar})

IT.insertXYZintoPLY(wallPoints_hmf,workingDirOF+'/Wall.ply',workingDirOF+'/WallDeformed.ply')
IT.insertXYZintoPLY(inletPoints_hmf,workingDirOF+'/Inlet.ply',workingDirOF+'/InletDeformed.ply')

plt.show()

