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
from copy import deepcopy

print 'Running HostMeshMotility...'

#Step one, create out contraction Matrix
workingDir = '/home/stephen/OpenFOAM/PhDSimulations/MotionBiPhasic/ContractionContours/'
workingDirOF = '/home/stephen/OpenFOAM/PhDSimulations/BiPhasicSimulation/'

#outfolder for hostmeshing details
os.mkdir(workingDirOF+'/HMresults')

#we will have a control file that has all of the contour groups we are planning to use
activeContours=open(workingDir+'ContourList').readlines()

#read in the contraction map
CMap=np.loadtxt(workingDir+'contractionMap')

#Set the number of processors
numProc=20

#Define the patch names
patchNames=['inlet','wall']
scale=0.001 #working in m in openFOAM, but HM is hardcoded for mm

# fititng parameters for host mesh fitting
host_mesh_pad = 20.0 # host mesh padding around slave points
host_elem_type = 'quad555' # quadrilateral cubic host elements
host_elems = [1,1,2] # a single element host mesh
maxit = 10
sobd = [15,15,15]
sobw = 2e-7 #0.000001#for quad444, 0.000005 is a good start, then crank
xtol = 1e-12

#===========================Static File IO==============================#
#set up a blank list that will store each group of patch nodes
#Loop through each processor, read the points and we assign them to a concated np array.
#we will store their index so we can seperate them laters
patchList=[0]*(numProc*len(patchNames))
print np.size(patchList)
for i in range(numProc):
	print '---Storing patches for processor',i,'---'
	for j in range(len(patchNames)):
		print j,((i*len(patchNames))+j)
		patchList[(i*len(patchNames))+j]=FT.readPatchFile(workingDirOF+'processor'+str(i)+'/constant/',patchNames[j]+'PatchNodes')

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
#passivePoints=passivePoints[::100]

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

#contractionPositions=open(workingDir+'MeshPositions').readlines()
#deformedMap=np.zeros(np.shape(contractionPositions)[0])
#for i in range(np.shape(contractionPositions)[0]):
#	deformedMap[i]=int(contractionPositions[i].rstrip().split()[1])
#	
#contractionMatrix,endtime=IT.createContractionMatrixDeformed(activeContours,contourInformation,contourData,CMap,deformedMap)

#lets make a nxm spatial array for the first time point, where n=(x,y,z) and m are the nodes
cMshapeIM=np.shape(contractionMatrix)
cMshape=[0]*4
cMshape[0]=np.shape(contractionMatrix)[0]
cMshape[1]=np.shape(contractionMatrix)[1]
cMshape[2]=np.shape(contractionMatrix[0][0])[0]
cMshape[3]=np.shape(contractionMatrix[0][0])[1]

#Create a list for passive nodes at each time position, to be used for splining
passiveList=np.zeros((passivePoints.shape[0],passivePoints.shape[1],cMshape[0]))
#======================================++++=================================#

#==========================Host Mesh Fitting================================#

#source points as contours
#temp=np.ones((cMshape[1],cMshape[2],cMshape[3]+1))
source_points_fitting=np.ones((1,3))
for i in range(cMshape[1]):
	cShape=np.shape(contractionMatrix[0][i])
	print 'Shape check:',cShape
	temp=np.ones((cShape[0],cShape[1]+1))
	temp[:,:2]=contractionMatrix[0][i]
	temp[:,2]=temp[:,2]*np.float(contourInformation[i][0])
	source_points_fitting=np.vstack([source_points_fitting,temp])

source_points_fitting=np.delete(source_points_fitting,(0),0)
print 'Source_points_fitting shape', np.shape(source_points_fitting)

#fixedContour1=np.loadtxt(workingDir+'FixedContours/C1')
#fixedContour2=np.loadtxt(workingDir+'FixedContours/C2')
#fixedContour3=np.loadtxt(workingDir+'FixedContours/C3')
#fixedContour4=np.loadtxt(workingDir+'FixedContours/C4')

#fixedContours=scipy.vstack([fixedContour1,fixedContour2,fixedContour3,fixedContour4])
#source_points_fitting=scipy.vstack([source_points_fitting,fixedContour1])


# make host mesh
host_mesh_orig = GFF.makeHostMeshMulti(
            	passivePoints.T,
                host_mesh_pad,
                host_elem_type,
                host_elems,
                )


#Store the original host mesh position
HostMeshOrig=host_mesh_orig.get_field_parameters()
test=host_mesh_orig.get_field_parameters()

print 'Host mesh shape'
print test.shape
hmfIndex=test.shape[1]/9
#test[2][3*hmfIndex:4*hmfIndex]=-400
#host_mesh_orig.set_field_parameters(test)
#print test[2][:,0]
test[2][1*hmfIndex:2*hmfIndex]=-670.5
test[2][2*hmfIndex:3*hmfIndex]=-607.0
test[2][3*hmfIndex:4*hmfIndex]=-543.5
test[2][4*hmfIndex:5*hmfIndex]=-480.0
test[2][5*hmfIndex:6*hmfIndex]=-451.75
test[2][6*hmfIndex:7*hmfIndex]=-423.5
test[2][7*hmfIndex:8*hmfIndex]=-395.25
#print test[2][:,0]
host_mesh_orig.set_field_parameters(test)

host_mesh=deepcopy(host_mesh_orig)


#embedding each patch For passive deform
passivePoints_xi = host_mesh.find_closest_material_points(
                            passivePoints,
                            initGD=[200,200,200],
                            verbose=True,
                            )[0]

passivePoints_passive = geometric_field.makeGeometricFieldEvaluatorSparse(
								host_mesh, [1,1],
    							matPoints=passivePoints_xi,
                                )

# calc slave node xi in host
print 'calculating slave xi...'
slave_xi = host_mesh.find_closest_material_points(
					source_points_fitting,
					initGD=[100,100,100],
					verbose=True,
					)[0]

#Loop through time points
for timePoint in range(cMshape[0]):

	print 'Processing time:', timePoint

	host_mesh=deepcopy(host_mesh_orig)

	#Read in the target points
	target_points=np.ones((1,3))
	for i in range(cMshape[1]):
		cShape=np.shape(contractionMatrix[0][i])
		print 'Shape check:',cShape
		temp=np.ones((cShape[0],cShape[1]+1))
		temp[:,:2]=contractionMatrix[timePoint][i]
		temp[:,2]=temp[:,2]*np.float(contourInformation[i][0])
		target_points=np.vstack([target_points,temp])

	target_points=np.delete(target_points,(0),0)
	print 'target_points shape', np.shape(target_points)
	#target_points=scipy.vstack([target_points,fixedContour1])

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
                	slave_xi=slave_xi,
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
	#we are also going to save host mesh data
	os.mkdir(workingDirOF+'/HMresults/'+str(timePoint+1))
	np.savetxt(workingDirOF+'/HMresults/'+str(timePoint+1)+'/target_points',target_points)
	np.savetxt(workingDirOF+'/HMresults/'+str(timePoint+1)+'/source_points_fitting',source_points_fitting)
	np.savetxt(workingDirOF+'/HMresults/'+str(timePoint+1)+'/source_points_fitting_hmf',source_points_fitting_hmf)
	np.savetxt(workingDirOF+'/HMresults/'+str(timePoint+1)+'/HostMeshOrig',HostMeshOrig[:,:,0].T)
	np.savetxt(workingDirOF+'/HMresults/'+str(timePoint+1)+'/HostMeshDeform',host_x_opt[:,:,0].T)

#=============================================================#
#Export the spine coefficients
passiveList=np.insert(passiveList,0,passivePoints,axis=2)
#Set up a time array (xn for splining)
timeVector=range(cMshape[0]+1)

pos=0

for p in range(numProc):
	os.mkdir(workingDirOF+'processor'+str(p)+'/constant/patchDisplacements')
	#os.mkdir(workingDirOF+'processor'+str(p)+'/constant/patchPositions')

	for i in range(len(patchNames)):
		print 'setting up file directories...'
		os.mkdir(workingDirOF+'processor'+str(p)+'/constant/patchDisplacements/'+str(patchNames[i]))
		#os.mkdir(workingDirOF+'processor'+str(p)+'/constant/patchPositions/'+str(patchNames[i]))
		for k in range(cMshape[0]+1):
			os.mkdir(workingDirOF+'processor'+str(p)+'/constant/patchDisplacements/'+str(patchNames[i])+'/'+str(k))
			#os.mkdir(workingDirOF+'processor'+str(p)+'/constant/patchPositions/'+str(patchNames[i])+'/'+str(k))

		if index[(p*len(patchNames))+i]!=-1:
			print 'creating spline files for'+patchNames[i]+'on proc'+str(p)
			patchNodes=passiveList[pos:index[(p*len(patchNames))+i],:,:]*scale
			#pos:index[(p*len(patchNames))+i]
		
			#So we will be making an a, b c and d vector files 
			aCoefs=np.zeros((np.shape(patchNodes)[0],3,cMshape[0]+1))
			bCoefs=np.zeros((np.shape(patchNodes)[0],3,cMshape[0]+1))
			cCoefs=np.zeros((np.shape(patchNodes)[0],3,cMshape[0]+1))
			dCoefs=np.zeros((np.shape(patchNodes)[0],3,cMshape[0]+1))
			#Loop through the number of nodes, and calculate the spines for x, y and z
			print 'Splining patch:', patchNames[i]
			for j in range(np.shape(patchNodes)[0]):
				for k in range(3): #number of dimensions
					yn=patchNodes[j,k,:]
					[a,b,c,d]=IT.basicCubicSpline(cMshape[0]+1,timeVector,yn)		
					[a,b,c,d]=IT.transformCubicCoeffiencets(a,b,c,d,timeVector)
			
					aCoefs[j,k,:]=a 
					bCoefs[j,k,:]=b
					cCoefs[j,k,:]=c
					dCoefs[j,k,:]=d
		
			#once we have filled the vectors, we write them to the time files
			print 'saving time files...', patchNames[i]
			#for time 0, we need our OF vector files
			#np.savetxt(workingDirOF+'processor'+str(p)+'/constant/patchPositions/'+str(patchNames[i])+'/'+str(0)+'/patchDisplacements',patchNodes[:,:,0])
			FT.writeVectorField(workingDirOF+'processor'+str(p)+'/constant/patchDisplacements/'+str(patchNames[i])+'/'+str(0)+'/'+str(patchNames[i])+'A',aCoefs[:,:,0])
			FT.writeVectorField(workingDirOF+'processor'+str(p)+'/constant/patchDisplacements/'+str(patchNames[i])+'/'+str(0)+'/'+str(patchNames[i])+'B',bCoefs[:,:,0])
			FT.writeVectorField(workingDirOF+'processor'+str(p)+'/constant/patchDisplacements/'+str(patchNames[i])+'/'+str(0)+'/'+str(patchNames[i])+'C',cCoefs[:,:,0])
			FT.writeVectorField(workingDirOF+'processor'+str(p)+'/constant/patchDisplacements/'+str(patchNames[i])+'/'+str(0)+'/'+str(patchNames[i])+'D',dCoefs[:,:,0])

			for k in range(1,cMshape[0]+1):
				np.savetxt(workingDirOF+'processor'+str(p)+'/constant/patchDisplacements/'+str(patchNames[i])+'/'+str(k)+'/a',aCoefs[:,:,k])
				np.savetxt(workingDirOF+'processor'+str(p)+'/constant/patchDisplacements/'+str(patchNames[i])+'/'+str(k)+'/b',bCoefs[:,:,k])
				np.savetxt(workingDirOF+'processor'+str(p)+'/constant/patchDisplacements/'+str(patchNames[i])+'/'+str(k)+'/c',cCoefs[:,:,k])
				np.savetxt(workingDirOF+'processor'+str(p)+'/constant/patchDisplacements/'+str(patchNames[i])+'/'+str(k)+'/d',dCoefs[:,:,k])
				#np.savetxt(workingDirOF+'processor'+str(p)+'/constant/patchPositions/'+str(patchNames[i])+'/'+str(k)+'/patchDisplacements',patchNodes[:,:,k])

			pos=index[(p*len(patchNames))+i]
		else:
			#just print a zero
			print 'printing a zero field'
			#FT.createBlankVectorField(workingDirOF+'processor'+str(p)+'/constant/patchPositions/'+str(patchNames[i])+'/'+str(0),'patchDisplacementsVector')
			FT.createBlankVectorField(workingDirOF+'processor'+str(p)+'/constant/patchDisplacements/'+str(patchNames[i])+'/'+str(0),str(patchNames[i])+'A')
			FT.createBlankVectorField(workingDirOF+'processor'+str(p)+'/constant/patchDisplacements/'+str(patchNames[i])+'/'+str(0),str(patchNames[i])+'B')
			FT.createBlankVectorField(workingDirOF+'processor'+str(p)+'/constant/patchDisplacements/'+str(patchNames[i])+'/'+str(0),str(patchNames[i])+'C')
			FT.createBlankVectorField(workingDirOF+'processor'+str(p)+'/constant/patchDisplacements/'+str(patchNames[i])+'/'+str(0),str(patchNames[i])+'D')

			for k in range(1,cMshape[0]+1):
				FT.spoofEmptyDisplacementFile(workingDirOF+'processor'+str(p)+'/constant/patchDisplacements/'+str(patchNames[i])+'/'+str(k)+'/a')
				FT.spoofEmptyDisplacementFile(workingDirOF+'processor'+str(p)+'/constant/patchDisplacements/'+str(patchNames[i])+'/'+str(k)+'/b')
				FT.spoofEmptyDisplacementFile(workingDirOF+'processor'+str(p)+'/constant/patchDisplacements/'+str(patchNames[i])+'/'+str(k)+'/c')
				FT.spoofEmptyDisplacementFile(workingDirOF+'processor'+str(p)+'/constant/patchDisplacements/'+str(patchNames[i])+'/'+str(k)+'/d')

	
#IT.basicCubicSplinePlot(a,b.tolist(),c.tolist(),d.tolist(),timeVector)
#plt.plot(timeVector,yn,'or')
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
