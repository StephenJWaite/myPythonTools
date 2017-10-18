#script reads in the motion files given to a OF simulation, and create a series of ply files for each given time intervals.
import FoamTools as FT
import numpy as np
import os,sys
import matplotlib.pyplot as plt

sys.path.insert(0,'./../')
workingDir='/home/stephen/OpenFOAM/Simulations2/Rumens/RumenSplineTest/'
patchName='inlet'
#Read in the mesh nodes
meshNodes=FT.readPatchFile(workingDir+'constant/',patchName+'PatchNodes')

#Set time range and dt
timeRange=5
dT=np.array([0.5,1])
numSteps=5

#Read in the spline coeffs for time 0
a=FT.readPatchFile(workingDir+'constant/patchDisplacements/'+patchName+'/0/',patchName+'A')
b=FT.readPatchFile(workingDir+'constant/patchDisplacements/'+patchName+'/0/',patchName+'B')
c=FT.readPatchFile(workingDir+'constant/patchDisplacements/'+patchName+'/0/',patchName+'C')
d=FT.readPatchFile(workingDir+'constant/patchDisplacements/'+patchName+'/0/',patchName+'D')

Nodes=np.zeros((meshNodes.shape[0],3))
dt=np.linspace(0,1,numSteps)
for i in range(timeRange-1):
	for t in range(dT.shape[0]):
		for node in range(meshNodes.shape[0]):
			Nodes[node,0] = a[node,0]*(dT[t]+i)**3 + b[node,0]*(dT[t]+i)**2 + c[node,0]*(dT[t]+i) + d[node,0]
			Nodes[node,1] = a[node,1]*(dT[t]+i)**3 + b[node,1]*(dT[t]+i)**2 + c[node,1]*(dT[t]+i) + d[node,1]
			Nodes[node,2] = a[node,2]*(dT[t]+i)**3 + b[node,2]*(dT[t]+i)**2 + c[node,2]*(dT[t]+i) + d[node,2]

		#now output the time file
		print Nodes.shape
		np.savetxt(workingDir+'/constant/SplineTest/'+patchName+'/'+str(i+dT[t])+'.xyz',Nodes)
	#Rea in the new abcd
	a=np.loadtxt(workingDir+'constant/patchDisplacements/'+patchName+'/'+str(i+1)+'/a')
	b=np.loadtxt(workingDir+'constant/patchDisplacements/'+patchName+'/'+str(i+1)+'/b')
	c=np.loadtxt(workingDir+'constant/patchDisplacements/'+patchName+'/'+str(i+1)+'/c')
	d=np.loadtxt(workingDir+'constant/patchDisplacements/'+patchName+'/'+str(i+1)+'/d')


#plt.plot(dt,x)
plt.show()


