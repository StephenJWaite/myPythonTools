#script is used to generate spline files for given mesh points.
import numpy as np
import imageTools as IT
workingDir='/home/stephen/OpenFOAM/Simulations2/Rumens/RumenSplineTest/constant/patchPositions/'
patchName='wall/'
outDir='/home/stephen/OpenFOAM/Simulations2/Rumens/RumenSplineTest/constant/manualSplineTests/'
meshDir='/home/stephen/OpenFOAM/Simulations2/Rumens/RumenSplineTest/constant/'

#read in the desired number of time steps
start=0
stop=1
template=np.loadtxt(workingDir+patchName+str(start)+'/patchDisplacements')
arraySize=[template.shape[0],template.shape[1],stop-start+1]
print arraySize
Displacements=np.zeros((arraySize[0],arraySize[1],arraySize[2]))
aC=np.zeros((arraySize[0],arraySize[1],arraySize[2]))
bC=np.zeros((arraySize[0],arraySize[1],arraySize[2]))
cC=np.zeros((arraySize[0],arraySize[1],arraySize[2]))
dC=np.zeros((arraySize[0],arraySize[1],arraySize[2]))

Displacements=np.zeros((arraySize[0],arraySize[1],arraySize[2]))

for i in range(start,stop+1):
	Displacements[:,:,i]=np.loadtxt(workingDir+patchName+str(i)+'/patchDisplacements')

#Now we each points, x y and z.
timeVector=range(start,stop+1)
for i in range(Displacements.shape[0]):
	for j in range(Displacements.shape[1]):
		[a,b,c,d]=IT.basicCubicSpline((stop-start+1),timeVector,Displacements[i,j,:])
		[a,b,c,d]=IT.transformCubicCoeffiencets(a,b,c,d,timeVector)
		aC[i,j,:]=a
		bC[i,j,:]=b
		cC[i,j,:]=c
		dC[i,j,:]=d

for i in range(len(timeVector)-1):
	np.savetxt(outDir+patchName+'a'+str(timeVector[i]),aC[:,:,i])
	np.savetxt(outDir+patchName+'b'+str(timeVector[i]),bC[:,:,i])
	np.savetxt(outDir+patchName+'c'+str(timeVector[i]),cC[:,:,i])
	np.savetxt(outDir+patchName+'d'+str(timeVector[i]),dC[:,:,i])



