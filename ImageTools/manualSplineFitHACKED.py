#script is used to generate spline files for given mesh points.
import numpy as np
import imageTools as IT
import FoamTools as FT
import os
workingDir='/home/stephen/OpenFOAM/Simulations2/Rumens/RumenSplineTest/constant/patchPositions/'
patchName='wall'
outDir='/home/stephen/OpenFOAM/Simulations2/Rumens/SplineChecks/RumenSplineTestCheck6/constant/patchDisplacements/'
meshDir='/home/stephen/OpenFOAM/Simulations2/Rumens/SplineChecks/RumenSplineTestCheck6/constant/'

#read in the desired number of time steps
start=0
stop=20
template=np.loadtxt(workingDir+patchName+'/'+str(start)+'/patchDisplacements')
arraySize=[template.shape[0],template.shape[1],stop-start+1]
print arraySize
Displacements=np.zeros((arraySize[0],arraySize[1],arraySize[2]))
aC=np.zeros((arraySize[0],arraySize[1],arraySize[2]))
bC=np.zeros((arraySize[0],arraySize[1],arraySize[2]))
cC=np.zeros((arraySize[0],arraySize[1],arraySize[2]))
dC=np.zeros((arraySize[0],arraySize[1],arraySize[2]))

#read in time 0
meshPos=FT.readPatchFile(meshDir,patchName+'PatchNodes') #hack line
Displacements=np.zeros((arraySize[0],arraySize[1],arraySize[2]))

Displacements[:,:,0]=meshPos
for i in range(start,stop):
	Displacements[:,:,i+1]=np.loadtxt(workingDir+patchName+'/'+str(i)+'/patchDisplacements')

#Now we each points, x y and z.
print 'Splinging, this will take a while'
timeVector=range(start,stop+1)
for i in range(Displacements.shape[0]):
	for j in range(Displacements.shape[1]):
		[a,b,c,d]=IT.basicCubicSpline((stop-start+1),timeVector,Displacements[i,j,:])
		[a,b,c,d]=IT.transformCubicCoeffiencets(a,b,c,d,timeVector)
		aC[i,j,:]=a
		bC[i,j,:]=b
		cC[i,j,:]=c
		dC[i,j,:]=d

print 'Writing out files'
os.mkdir(outDir+patchName)
for i in range(len(timeVector)-1):
	os.mkdir(outDir+patchName+'/'+str(timeVector[i]))
	np.savetxt(outDir+patchName+'/'+str(timeVector[i])+'/a',aC[:,:,i])
	np.savetxt(outDir+patchName+'/'+str(timeVector[i])+'/b',bC[:,:,i])
	np.savetxt(outDir+patchName+'/'+str(timeVector[i])+'/c',cC[:,:,i])
	np.savetxt(outDir+patchName+'/'+str(timeVector[i])+'/d',dC[:,:,i])

#write out those mesh positions plox
#np.savetxt(outDir+patchName+'/patchDisplacement/0/patchDisplacement',Displacements[:,:,0])
#np.savetxt(outDir+patchName+'/patchDisplacement/1/patchDisplacement',Displacements[:,:,1])



