#visualise patch points for sanity
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

plt.switch_backend('TkAgg')

workingDir='/home/stephen/OpenFOAM/Simulations2/Rumens/NewRumen/constant'
wallPatch=np.loadtxt(workingDir+'/wallPatchNodes_outfile.xyz')
inletPatch=np.loadtxt(workingDir+'/inletPatchNodes_outfile.xyz')
outletPatch=np.loadtxt(workingDir+'/outletPatchNodes_outfile.xyz')

fig=plt.figure(2)
visPlot=fig.add_subplot(111,projection='3d')
visPlot.set_xlim(50,350)
visPlot.set_ylim(100,400)
visPlot.set_zlim(-700,-400)
visPlot.plot(inletPatch[:,0],inletPatch[:,1],inletPatch[:,2],'ob')
visPlot.plot(wallPatch[:,0],wallPatch[:,1],wallPatch[:,2],'xr')

plt.show()