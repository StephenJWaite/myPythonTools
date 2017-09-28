#This is a basic host mesh fitting script to test componants.
import sys
sys.path.insert(0,'./../../lib')
import numpy as np
import scipy
import itertools
import copy
import imageTools as IT
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import cKDTree
from fieldwork.field.tools import fitting_tools
from fieldwork.field import geometric_field
from fieldwork.field import geometric_field_fitter as GFF
from gias.common import fieldvi, transform3D
from gias.common import alignment_fitting as af


print 'Running HostMeshMotility...'

#Step one, create out contraction Matrix
workingDir = './../imageToolsTestSpace/'

#we will have a control file that has all of the contour groups we are planning to use
activeContours=open(workingDir+'ContourList').readlines()

#read in the contraction map
CMap=np.loadtxt(workingDir+'contractionMap')

#set up infromation arrays
contourData=[]
contourInformation=[]
#loop through activeContours
for line in activeContours:
	print 'Opening contour:',line.rstrip()
	zPos,numContours,contours=IT.readContourFile(workingDir+line.rstrip()+'/Results.txt')
	#store the zPos and numContous in contourInformation
	contourInformation=contourInformation+[[zPos,numContours]]
	contourData=contourData+[contours]

contractionMatrix,endtime=IT.createContractionMatrix(activeContours,contourInformation,contourData,CMap)

#lets make a nxm spatial array for the first time point, where n=(x,y,z) and m are the nodes
cMshape=np.shape(contractionMatrix)

#Set up plot variables
fig2=plt.figure(2)
timePlot=fig2.add_subplot(111,projection='3d')
timePlot.set_xlim(50,350)
timePlot.set_ylim(100,400)
plt.ion()

#=============================================================================#
# fititng parameters for host mesh fitting
host_mesh_pad = 10.0 # host mesh padding around slave points
host_elem_type = 'quad333' # quadrilateral cubic host elements
host_elems = [1,1,3] # a single element host mesh
maxit = 10
sobd = [4,4,4]
sobw = 0.00001#1e-10
xtol = 1e-12

# source points for fitting
source_points_fitting_file = workingDir+'RumenWithInletOutlet_outfile.xyz'
source_points_passive = np.loadtxt(source_points_fitting_file)

#source points as contours
temp=np.ones((cMshape[1],cMshape[2],cMshape[3]+1))
for i in range(cMshape[1]):
	print 'Shape check:',np.shape(temp[i,:,:2]),np.shape(contractionMatrix[0][i])
	temp[i,:,:2]=(contractionMatrix[0][i])-256
	print 'zPos:',contourInformation[i][0]
	temp[i,:,2]=temp[i,:,2]*np.float(contourInformation[i][0])
	#temp[i,:,2]=temp[i,:,2]*-np.float(contourInformation[i][0])#-375

#Read in any fixed points
fixedContour1=np.loadtxt(workingDir+'FixedContours/InletFixedContour1')
fixedContour2=np.loadtxt(workingDir+'FixedContours/InletFixedContour2')
fixedContour3=np.loadtxt(workingDir+'FixedContours/FixedContour1')
fixedContour4=np.loadtxt(workingDir+'FixedContours/OutletFixedContour1')
fixedContour5=np.loadtxt(workingDir+'FixedContours/OutletFixedContour2')

fixedContours=scipy.vstack([fixedContour1,fixedContour2,fixedContour3,fixedContour4,fixedContour5])

print 'countour shape:',np.shape(temp.reshape((cMshape[1]*cMshape[2],cMshape[3]+1)))
#print 'fixedcontour shape',fixedContour.shape

source_points_fitting = temp.reshape((cMshape[1]*cMshape[2],cMshape[3]+1))

#With fixed contours
#source_points_fitting = scipy.vstack([temp.reshape((cMshape[1]*cMshape[2],cMshape[3]+1)),fixedContours])

#########################################################################################################
#Loop through time points
for timePoint in range(11,12):#range(20,cMshape[0]):

	print 'fitting time: ',timePoint
	temp=np.ones((cMshape[1],cMshape[2],cMshape[3]+1))
	for i in range(cMshape[1]):
		print 'Shape check:',np.shape(temp[i,:,:2]),np.shape(contractionMatrix[timePoint][i])
		temp[i,:,:2]=(contractionMatrix[timePoint][i])-256
		print 'zPos:',contourInformation[i][0]
		temp[i,:,2]=temp[i,:,2]*np.float(contourInformation[i][0])
		#temp[i,:,2]=temp[i,:,2]*-np.float(contourInformation[i][0])-375

	timePlot.clear()
	timePlot.set_xlim(50,350)
	timePlot.set_ylim(100,400)
	timePlot.view_init(elev=-81+i, azim=-113+i)
	timePlot.plot(contractionMatrix[0][0][:,0],contractionMatrix[0][0][:,1],5,'-xb')
	timePlot.plot(contractionMatrix[0][1][:,0],contractionMatrix[0][1][:,1],1,'-xb')
	timePlot.plot(contractionMatrix[timePoint][0][:,0],contractionMatrix[timePoint][0][:,1],5,'-or')
	timePlot.plot(contractionMatrix[timePoint][1][:,0],contractionMatrix[timePoint][1][:,1],1,'-or')
	plt.pause(0.05)
	# source points to be passived deformed (not fitted)
	#target_points_file = '/home/stephen/Documents/PhD_data/RandomSTLS/RumenSlices2.xyz'
	target_points = temp.reshape((cMshape[1]*cMshape[2],cMshape[3]+1))

	#With fixed contours
	#target_points = scipy.vstack([temp.reshape((cMshape[1]*cMshape[2],cMshape[3]+1)),fixedContours])

	# host mesh fit source fitting points to target points and
	# apply HMF transform to passive source points

	# define some slave obj funcs
	target_tree = cKDTree(target_points)

	# distance between each target point and its closest source fitting point
	# should not use if source has less geometry than target
	def slave_func_tpsp(x):
		#print 'Slave_Func check:x.shape ',np.shape(x)
		sourcetree = cKDTree(x)
		d = sourcetree.query(target_points)[0]
		#print 'Slave_Func check:d ',np.shape(d) <- this has this shape, because its calcuating the closest source points for each target point
		return d

	print 'source shape: ',np.shape(source_points_fitting)
	print 'target shape: ',np.shape(target_points)
	#my super disgusting self coded distance function cause i couldnt find a scipy or numpy tool to do it...
	#also its expecting the shape to be points,coords
	def my_slave_func(v2):
		v1=target_points
		d=np.sqrt((v1[:,0]-v2[:,0])**2 + (v1[:,1]-v2[:,1])**2 + (v1[:,2]-v2[:,2])**2)
		return d

	slave_func = my_slave_func

	# make host mesh
	host_mesh = GFF.makeHostMeshMulti(
            		source_points_passive.T,
                	host_mesh_pad,
                	host_elem_type,
                	host_elems,
                	)


	#Store the original host mesh position
	HostMeshOrig=host_mesh.get_field_parameters()

	#For passive deform
	# calculate the emdedding (xi) coordinates of passive
	# source points.
	source_points_passive_xi = host_mesh.find_closest_material_points(
                            	source_points_passive,
                            	initGD=[50,50,50],
                            	verbose=True,
                            	)[0]

	#make passive source point evaluator function
	eval_source_points_passive = geometric_field.makeGeometricFieldEvaluatorSparse(
									host_mesh, [1,1],
    								matPoints=source_points_passive_xi,
                                	)

	# host mesh fit
	host_x_opt, source_points_fitting_hmf,\
	slave_xi, rmse_hmf = fitting_tools.hostMeshFitPoints(
                	host_mesh,
                	source_points_fitting,
                	slave_func,
                	max_it=maxit,
                	sob_d=sobd,
                	sob_w=sobw,
                	verbose=True,
                	xtol=xtol
                    )
	
	# evaluate the new positions of the passive source points
	source_points_passive_hmf = eval_source_points_passive(host_x_opt).T

	#write out the Results
	np.savetxt(workingDir+'./../Geoms/SoblovTests/SB'+str(timePoint)+'.xyz',source_points_passive_hmf)
#########################################################################################################



#plt.ioff()
#plt.show()
# view
v = fieldvi.fieldvi()
#v.addData('Origin',np.zeros((4,3)),renderArgs={'mode':'sphere', 'scale_factor':5, 'color':(1,0,0)})
v.addData('target points', target_points, renderArgs={'mode':'sphere', 'scale_factor':3, 'color':(1,0,0)})
v.addData('source points fitting points', source_points_fitting, renderArgs={'mode':'sphere', 'scale_factor':2, 'color':(0.5,0,0)})
v.addData('source points fitting', source_points_fitting, renderArgs={'mode':'sphere', 'scale_factor':3, 'color':(3,0,0)})
v.addData('Host Mesh Orig', HostMeshOrig[:,:,0].T, renderArgs={'mode':'sphere','scale_factor':5, 'color':(0,0,0.5)})
v.addData('source points fitting hmf', source_points_fitting_hmf, renderArgs={'mode':'point'})
v.addData('Host Mesh Deformed', host_x_opt[:,:,0].T, renderArgs={'mode':'sphere','scale_factor':5, 'color':(0,0.25,0)})
v.addData('source points passive', source_points_passive, renderArgs={'mode':'sphere', 'scale_factor':1, 'color':(0,0,0.6)})
v.addData('source points passive hmf', source_points_passive_hmf, renderArgs={'mode':'sphere', 'scale_factor':1, 'color':(0,0,0.6)})

v.configure_traits()
v.scene.background=(0,0,0)
