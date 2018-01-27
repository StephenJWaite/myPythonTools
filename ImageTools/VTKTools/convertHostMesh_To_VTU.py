#convertHostMesh_To_VTU.py converts a host mesh node file into 
#a vtu file wit connectivity, word
from pyevtk.hl import unstructuredGridToVTK
from pyevtk.vtk import VtkTriangle, VtkQuad
import numpy as np

#read in the point file.
points=np.loadtxt('hostMeshPoints')
print np.shape(points)
print points

#deine the nx ny nz
nx=1
ny=1
nz=2
#elements are stored
#  _____
# |1_|2_|
# |3_|4_|
#1   2  3
#

#construct elements
elems=np.zeros((nx*ny*nz,8))
stencil=np.zeros(8)
stencil=np.array([1,2,4,5,10,11,13,14])
stencil=np.array([0,1,3,4,9,10,12,13])
elem2=[]

for i in range(0,nx):
	for j in range(0,ny):
		for k in range(0,nz):
			#print stencil+(i+(nx+1)*j)+((nx+1)*(ny+1))*k
			print i,j,k
			#print k+(ny+1)*j+(ny+1)*(nx+1)*i #holy shit this is the stencil
			print k+ny*j+ny*nx*i
			elems[(k+ny*j+ny*nx*i)]=stencil+(i+(nx+1)*j)+((nx+1)*(ny+1))*k

print elems

#conn = np.zeros(nx*ny*nz*2*4)

x=np.asfortranarray(points[:,0])
y=np.asfortranarray(points[:,1])
z=np.asfortranarray(points[:,2])

conn=np.zeros(8)
conn[0],conn[1],conn[2],conn[3]=0,1,4,3
conn[4],conn[5],conn[6],conn[7]=9,10,13,12

offset = np.zeros(2)
offset[0] = 4
offset[1] = 8


ctype = np.zeros(2)
ctype[0], ctype[1] = VtkQuad.tid, VtkQuad.tid



#unstructuredGridToVTK("moose", x, y, z, connectivity = conn, offsets = offset, cell_types = ctype, cellData = None, pointData = None)
elemOffset=[]
elemCon=[]
elemcType=[]
#Construct and elementt, they are hexes
for i in range(0,elems.shape[0]):
	#face 1: 1,3,4,2 front
	elemCon.append(elems[i,(0,2,3,1)])
	elemOffset.append((i*24)+4)
	elemcType.append(VtkQuad.tid)
	#face 2: 2,4,8,6 right
	elemCon.append(elems[i,(1,3,7,5)])	
	elemOffset.append((i*24)+8)
	elemcType.append(VtkQuad.tid)
	#face 3: 1,3,7,5 left
	elemCon.append(elems[i,(0,2,6,4)])
	elemOffset.append((i*24)+12)
	elemcType.append(VtkQuad.tid)
	#face 4: 5,7,8,6 back
	elemCon.append(elems[i,(4,6,7,5)])
	elemOffset.append((i*24)+16)
	elemcType.append(VtkQuad.tid)
	#face 5: 3,7,8,4 top
	elemCon.append(elems[i,(2,6,7,3)])
	elemOffset.append((i*24)+20)
	elemcType.append(VtkQuad.tid)
	#face 6: 1,5,6,2 bottom
	elemCon.append(elems[i,(0,4,5,1)])
	elemOffset.append((i*24)+24)
	elemcType.append(VtkQuad.tid)


dims=np.shape(elemCon)
elemCon=np.asarray(elemCon).reshape(dims[0]*dims[1])
elemOffset=np.asarray(elemOffset)
elemcType=np.asarray(elemcType)

print 'elemCon', np.shape(elemCon)
print elemCon
print elemCon.size
print 'elemcType'
print elemcType
print elemcType.size
print 'elemOffset'
print elemOffset
print elemOffset.size

unstructuredGridToVTK("moose2", x, y, z, connectivity = elemCon, offsets = elemOffset, cell_types = elemcType, cellData = None, pointData = None)

