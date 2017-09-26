#Scrip is used to apply hostMesh fitting to parallised OF files, for any amount of patchs. Dope.
import os,sys
sys.path.insert(0,'./../../lib')
import numpy as np
import scipy
import itertools
import copy
import imageTools as IT
import FoamTools as FT
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import cKDTree
from fieldwork.field.tools import fitting_tools
from fieldwork.field import geometric_field
from fieldwork.field import geometric_field_fitter as GFF
from gias.common import fieldvi, transform3D
from gias.common import alignment_fitting as af

print 'Running HostMeshMotility...'

#Set the working directory
workingDirOF = '/home/stephen/OpenFOAM/Simulations2/Rumens/RumenParallelTest/'

#Importing the FOAM files as passive slave points
#Set the number of processors
numProc=4
patchNames=['inlet','outlet','wall']
scale=0.001 #working in m now not mm

#Not, this is hardcoded to work with three patches, wall, inlet and outlet, but its trivial to 
#extend it to any patch grou you want

#set up a blank list that will store each group of patch nodes
#Loop through each processor, read the points and we assign them to a concated np array.
#we will store their index so we can seperate them laters
patchList=[0]*(numProc*3)
print np.size(patchList)
for i in range(numProc):
	print '---Storing patches for processor',i,'---'
	for j in range(len(patchNames)):
		print j,((i*(numProc-1))+j)
		patchList[(i*(numProc-1))+j]=FT.readPatchFile(workingDirOF+'processor'+str(i)+'/constant/',patchNames[j]+'PatchNodes')

#Now loop through all the patchList, and smoosh it all into one np array, and store the start point of each patch.
index=[]
currentIndex=0
print np.shape(patchList)
passivePoints=np.zeros((1,3))
for patch in patchList:
	print np.shape(patch)
	#check if the patch isnt in a given processor.
	if len(patch)==0:
		print 'empty patch'
		index=index+[-1]
	else:
		currentIndex=currentIndex+np.shape(patch)[0]
		index=index+[currentIndex]
		passivePoints=np.vstack([passivePoints,np.asarray(patch)])

		
print index

passivePoints=np.delete(passivePoints,(0),0)/scale


#=============================================================#
#Now we can host mesh fit the points as normal
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
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

# fititng parameters for host mesh fitting
host_mesh_pad = 10.0 # host mesh padding around slave points
host_elem_type = 'quad333' # quadrilateral cubic host elements
host_elems = [1,1,2] # a single element host mesh
maxit = 10
sobd = [4,4,4]
sobw = 1e-10
xtol = 1e-12

#source points as contours
temp=np.ones((cMshape[1],cMshape[2],cMshape[3]+1))
for i in range(cMshape[1]):
	print 'Shape check:',np.shape(temp[i,:,:2]),np.shape(contractionMatrix[0][i])
	temp[i,:,:2]=(contractionMatrix[0][i])-256
	print 'zPos:',contourInformation[i][0]
	temp[i,:,2]=temp[i,:,2]*-np.float(contourInformation[i][0])-375

print 'countour shape:',np.shape(temp.reshape((cMshape[1]*cMshape[2],cMshape[3]+1)))

source_points_fitting = temp.reshape((cMshape[1]*cMshape[2],cMshape[3]+1))

#Loop through time points
for timePoint in range(11,12):#range(20,cMshape[0]):

	print 'fitting time: ',timePoint
	temp=np.ones((cMshape[1],cMshape[2],cMshape[3]+1))
	for i in range(cMshape[1]):
		print 'Shape check:',np.shape(temp[i,:,:2]),np.shape(contractionMatrix[timePoint][i])
		temp[i,:,:2]=(contractionMatrix[timePoint][i])-256
		print 'zPos:',contourInformation[i][0]
		temp[i,:,2]=temp[i,:,2]*-np.float(contourInformation[i][0])-375

	#Set up the target points
	target_points = temp.reshape((cMshape[1]*cMshape[2],cMshape[3]+1))

	# define some slave obj funcs
	target_tree = cKDTree(target_points)

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
            		passivePoints.T,
                	host_mesh_pad,
                	host_elem_type,
                	host_elems,
                	)


	#Store the original host mesh position
	HostMeshOrig=host_mesh.get_field_parameters()

	#embedding each patch For passive deform
	passivePoints_xi = host_mesh.find_closest_material_points(
                            	passivePoints,
                            	initGD=[50,50,50],
                            	verbose=True,
                            	)[0]

	passivePoints_passive = geometric_field.makeGeometricFieldEvaluatorSparse(
									host_mesh, [1,1],
    								matPoints=passivePoints_xi,
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
	passivePoints_hmf = passivePoints_passive(host_x_opt).T


#Now we are going to break up and write out these points back into patches.
#passivePoints_hmf=passivePoints
newPatch=[0]*(numProc*3)
pos=0
for i in range(numProc):
	os.mkdir(workingDirOF+'processor'+str(i)+'/constant/patchDisplacements')
	for j in range(len(patchNames)):
		if index[(i*(numProc-1))+j]!=-1:
			#we print the patchArrray at these points
			print index[(i*(numProc-1))+j]
			print pos,index[(i*(numProc-1))+j]
			np.savetxt(workingDirOF+'processor'+str(i)+'/constant/patchDisplacements/'+patchNames[j]+'Displacement',passivePoints_hmf[pos:index[(i*(numProc-1))+j],:]*scale)
			FT.writeXYZtoPointsVectorField(workingDirOF+'processor'+str(i)+'/constant/patchDisplacements/',patchNames[j])
			pos=index[(i*(numProc-1))+j]
		else:
			#just print a zero
			print i,j,patchNames[j]
			print [workingDirOF+'processor'+str(i)+'/constant/patchDisplacements/'+patchNames[j]+'Displacement']
			fid=open(workingDirOF+'processor'+str(i)+'/constant/patchDisplacements/'+patchNames[j]+'Displacement','w')
			fid.write('{0}'.format(0))
			fid.close()
			FT.createBlankVectorField(workingDirOF+'processor'+str(i)+'/constant/patchDisplacements/',patchNames[j])

# view
v = fieldvi.fieldvi()
#v.addData('Origin',np.zeros((4,3)),renderArgs={'mode':'sphere', 'scale_factor':5, 'color':(1,0,0)})
v.addData('target points', target_points, renderArgs={'mode':'sphere', 'scale_factor':3, 'color':(1,0,0)})
v.addData('source points fitting points', source_points_fitting, renderArgs={'mode':'sphere', 'scale_factor':2, 'color':(0.5,0,0)})
#v.addData('source points fitting', source_points_fitting, renderArgs={'mode':'sphere', 'scale_factor':3, 'color':(3,0,0)})
v.addData('Host Mesh Orig', HostMeshOrig[:,:,0].T, renderArgs={'mode':'sphere','scale_factor':5, 'color':(0,0,0.5)})
v.addData('source points fitting hmf', source_points_fitting_hmf, renderArgs={'mode':'point'})
v.addData('Host Mesh Deformed', host_x_opt[:,:,0].T, renderArgs={'mode':'sphere','scale_factor':5, 'color':(0,0.25,0)})
v.addData('patch points hmf',passivePoints_hmf, renderArgs={'mode':'sphere', 'scale_factor':1, 'color':(0,0,0.2)})


v.configure_traits()
v.scene.background=(0,0,0)

	

	

	

