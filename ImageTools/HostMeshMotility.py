#This is a basic host mesh fitting script to test componants.
import sys
sys.path.insert(0,'/hpc/swai013/Python/lib')
import numpy as np
import itertools
import copy
import imageTools as IT
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
temp=np.ones((cMshape[1],cMshape[2],cMshape[3]+1))
for i in range(cMshape[1]):
	print 'Shape check:',np.shape(temp[i,:,:2]),np.shape(contractionMatrix[0][i])
	temp[i,:,:2]=(contractionMatrix[0][i])-256
	print 'zPos:',contourInformation[i][0]
	temp[i,:,2]=temp[i,:,2]*-np.float(contourInformation[i][0])-375

#=============================================================================#
# fititng parameters for host mesh fitting
host_mesh_pad = 1.0 # host mesh padding around slave points
host_elem_type = 'quad333' # quadrilateral cubic host elements
host_elems = [1,1,2] # a single element host mesh
maxit = 10
sobd = [4,4,4]
sobw = 1e-10
xtol = 1e-12
 
# source points for fitting
source_points_fitting_file = workingDir+'Test.xyz'
source_points_fitting = np.loadtxt(source_points_fitting_file)
 
# source points to be passived deformed (not fitted)
#target_points_file = '/home/stephen/Documents/PhD_data/RandomSTLS/RumenSlices2.xyz'
target_points = temp.reshape((cMshape[1]*cMshape[2],cMshape[3]+1))

# host mesh fit source fitting points to target points and
# apply HMF transform to passive source points
 
# define some slave obj funcs
target_tree = cKDTree(target_points)

# distance between each target point and its closest source fitting point
# should not use if source has less geometry than target
def slave_func_tpsp(x):
    sourcetree = cKDTree(x)
    d = sourcetree.query(target_points)[0]
    return d
 
slave_func = slave_func_tpsp

# make host mesh
host_mesh = GFF.makeHostMeshMulti(
                source_points_fitting.T,
                host_mesh_pad,
                host_elem_type,
                host_elems,
                )
 
 
#Store the original host mesh position
HostMeshOrig=host_mesh.get_field_parameters()

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

#=============================================================#
# view
v = fieldvi.fieldvi()
v.addData('Origin',np.zeros((4,3)),renderArgs={'mode':'sphere', 'scale_factor':5, 'color':(1,0,0)})
v.addData('target points', target_points, renderArgs={'mode':'sphere', 'scale_factor':5, 'color':(1,0,0)})
v.addData('source points fitting', source_points_fitting, renderArgs={'mode':'point'})
v.addData('Host Mesh Orig', HostMeshOrig[:,:,0].T, renderArgs={'mode':'sphere','scale_factor':5, 'color':(0,0,0.5)})
v.addData('source points fitting hmf', source_points_fitting_hmf, renderArgs={'mode':'point'})
v.addData('Host Mesh Deformed', host_x_opt[:,:,0].T, renderArgs={'mode':'sphere','scale_factor':1.25, 'color':(0,0.25,0)})

v.configure_traits()
v.scene.background=(0,0,0)