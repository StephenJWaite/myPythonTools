#script reads in the motion files given to a OF simulation, and create a series of ply files for each given time intervals.
import FoamTools as FT
import numpy as np
import os,sys
sys.path.insert(0,'./../')
workingDir='/home/stephen/OpenFOAM/Simulations2/Rumens/RumenSplineTest/'

#Read in the mesh nodes
meshNodes=FT.readPatchFile(workingDir+'constant/','inletPatchNodes')

#Set time range and dt
timeRange=5
numSteps=5

#Read in the spline coeffs for time 0
a=FT.readPatchFile(workingDir+'constant/patchDisplacements/inlet/0/','inletA')
b=FT.readPatchFile(workingDir+'constant/patchDisplacements/inlet/0/','inletB')
c=FT.readPatchFile(workingDir+'constant/patchDisplacements/inlet/0/','inletC')
d=FT.readPatchFile(workingDir+'constant/patchDisplacements/inlet/0/','inletD')

x=np.zeros((5))
dt=np.linspace(0,1,numTime)
for t in range(numSteps):
	x[t] = a[0,0]*dt[t]**3 + a[0,0]*dt[t]**2 + a[0,0]*dt[t] + d[0,0]

print plotting
plt.plot(dt,x)
plt.show()


