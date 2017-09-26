#scaleXYZ is a script to scale hostmesh fit files oops
import numpy as np
workingDir='/home/stephen/OpenFOAM/Simulations2/Rumens/NewRumenNesiSetUp/constant/patchDisplacements/'
fileName='wall'
#Read in xyz file
point=np.loadtxt(workingDir+'unScaled/'+fileName+'Displacement')
#scale factor 0.001
scale=0.001
transformedPoints=point*scale
#write the file out
np.savetxt(workingDir+fileName+'Displacement',transformedPoints)

