#This is a basic host mesh fitting script to test componants.
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

#Step one, create out contraction Matrix
workingDir = './../imageToolsTestSpace/'
workingDirOF = '/home/stephen/OpenFOAM/Simulations2/Rumens/RumenBMES'

#we will have a control file that has all of the contour groups we are planning to use
activeContours=open(workingDir+'ContourList').readlines()

#read in the contraction map
CMap=np.loadtxt(workingDir+'contractionMap')

#Define the patch names
patchNames=['inlet','wall']
scale=0.001 #working in m in openFOAM, but HM is hardcoded for mm

# fititng parameters for host mesh fitting
host_mesh_pad = 10.0 # host mesh padding around slave points
host_elem_type = 'quad333' # quadrilateral cubic host elements
host_elems = [1,1,1] # a single element host mesh
maxit = 10
sobd = [4,4,4]
sobw = 0.00001#1e-10
xtol = 1e-12

#===========================Static File IO==============================#
#Read in patch files, and form the passive node array.
#load in the OF patches
patchList=[0]*len(patchNames)
for i in range(len(patchNames)):
	patchList[i]=FT.readPatchFile(workingDirOF+'/constant/',patchNames[i]+'PatchNodes')

index=[]
currentIndex=0
print np.shape(patchList)
passivePoints=np.zeros((1,3))
for patch in patchList:
	currentIndex=currentIndex+np.shape(patch)[0]
	index=index+[currentIndex]
	passivePoints=np.vstack([passivePoints,np.asarray(patch)])[0::5]

print index

passivePoints=np.delete(passivePoints,(0),0)/scale

#set up contour infromation arrays
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
cMshapeIM=np.shape(contractionMatrix)
cMshape=[0]*4
cMshape[0]=4
cMshape[1]=cMshapeIM[1]
cMshape[2]=cMshapeIM[2]
cMshape[3]=cMshapeIM[3]

#Create a list for passive nodes at each time position, to be used for splining
passiveList=np.zeros((passivePoints.shape[0],passivePoints.shape[1],cMshape[0]))
#======================================++++=================================#

#==========================Host Mesh Fitting================================#

#source points as contours
temp=np.ones((cMshape[1],cMshape[2],cMshape[3]+1))
for i in range(cMshape[1]):
	print 'Shape check:',np.shape(temp[i,:,:2]),np.shape(contractionMatrix[0][i])
	temp[i,:,:2]=(contractionMatrix[0][i])-256
	print 'zPos:',contourInformation[i][0]
	temp[i,:,2]=temp[i,:,2]*np.float(contourInformation[i][0])

print 'countour shape:',np.shape(temp.reshape((cMshape[1]*cMshape[2],cMshape[3]+1)))

source_points_fitting = temp.reshape((cMshape[1]*cMshape[2],cMshape[3]+1))



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

#Loop through time points
for timePoint in range(cMshape[0]):

	#Read in the target points
	print 'fitting time: ',timePoint
	temp=np.ones((cMshape[1],cMshape[2],cMshape[3]+1))
	for i in range(cMshape[1]):
		print 'Shape check:',np.shape(temp[i,:,:2]),np.shape(contractionMatrix[timePoint][i])
		temp[i,:,:2]=(contractionMatrix[timePoint][i])-256
		print 'zPos:',contourInformation[i][0]
		temp[i,:,2]=temp[i,:,2]*np.float(contourInformation[i][0])

	# source points to be passived deformed (not fitted)
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

	#put the passivePoints into a list
	passiveList[:,:,timePoint]=passivePoints_hmf

#=============================================================#
#Export the spine coefficients
#Set up a time array (xn for splining)
timeVector=range(cMshape[0])

pos=0
os.mkdir(workingDirOF+'/constant/patchDisplacements')
for i in range(len(patchNames)):
	os.mkdir(workingDirOF+'/constant/patchDisplacements/'+str(patchNames[i]))
	patchNodes=passiveList[pos:index[i],:,:]*scale
	#So we will be making an a, b c and d vector files 
	aCoefs=np.zeros((np.shape(patchNodes)[0],3,cMshape[0]))
	bCoefs=np.zeros((np.shape(patchNodes)[0],3,cMshape[0]))
	cCoefs=np.zeros((np.shape(patchNodes)[0],3,cMshape[0]))
	dCoefs=np.zeros((np.shape(patchNodes)[0],3,cMshape[0]))
	#Loop through the number of nodes, and calculate the spines for x, y and z
	for j in range(np.shape(patchNodes)[0]):
		for k in range(3): #number of dimensions
			yn=patchNodes[j,k,:]
			[a,b,c,d]=IT.basicCubicSpline(cMshape[0],timeVector,yn)
			[a,b,c,d]=IT.transformCubicCoeffiencets(a,b,c,d,timeVector)
			aCoefs[j,k,:]=a 
			bCoefs[j,k,:]=b
			cCoefs[j,k,:]=c
			dCoefs[j,k,:]=d
		
	#ones we have filled the vectors, we write them to the time files
	for k in range(cMshape[0]):
		os.mkdir(workingDirOF+'/constant/patchDisplacements/'+str(patchNames[i])+'/'+str(k))
		np.savetxt(workingDirOF+'/constant/patchDisplacements/'+str(patchNames[i])+'/'+str(k)+'/a',aCoefs[:,:,k])
		np.savetxt(workingDirOF+'/constant/patchDisplacements/'+str(patchNames[i])+'/'+str(k)+'/b',bCoefs[:,:,k])
		np.savetxt(workingDirOF+'/constant/patchDisplacements/'+str(patchNames[i])+'/'+str(k)+'/c',cCoefs[:,:,k])
		np.savetxt(workingDirOF+'/constant/patchDisplacements/'+str(patchNames[i])+'/'+str(k)+'/d',dCoefs[:,:,k])
	
	#for j in range(cMshape[0]):
	#	np.savetxt(workingDirOF+'/constant/patchDisplacements/'+patchNames[i]+str(timeVector[j]),patchNodes[:,:,j])

	pos=index[i]

	
IT.basicCubicSplinePlot(a,b.tolist(),c.tolist(),d.tolist(),timeVector)
plt.plot(timeVector,yn,'or')
plt.show()



# view
#v = fieldvi.fieldvi()
#v.addData('Origin',np.zeros((4,3)),renderArgs={'mode':'sphere', 'scale_factor':5, 'color':(1,0,0)})
#v.addData('target points', target_points, renderArgs={'mode':'sphere', 'scale_factor':3, 'color':(1,0,0)})
#v.addData('source points fitting points', source_points_fitting, renderArgs={'mode':'sphere', 'scale_factor':2, 'color':(0.5,0,0)})
#v.addData('source points fitting', source_points_fitting, renderArgs={'mode':'sphere', 'scale_factor':3, 'color':(3,0,0)})
#v.addData('Host Mesh Orig', HostMeshOrig[:,:,0].T, renderArgs={'mode':'sphere','scale_factor':5, 'color':(0,0,0.5)})
#v.addData('source points fitting hmf', source_points_fitting_hmf, renderArgs={'mode':'point'})
#v.addData('Host Mesh Deformed', host_x_opt[:,:,0].T, renderArgs={'mode':'sphere','scale_factor':5, 'color':(0,0.25,0)})
#v.addData('passive points', passivePoints, renderArgs={'mode':'sphere', 'scale_factor':1, 'color':(0,0,0.2)})
#v.addData('passive points hmf', passivePoints_hmf, renderArgs={'mode':'sphere', 'scale_factor':1, 'color':(0,0,0.5)})


#v.configure_traits()
#v.scene.background=(0,0,0)
