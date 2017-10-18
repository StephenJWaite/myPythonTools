#script is used to compare the linear equation for splines to 2 mesh positions, a quick sanity check
#before I check this in OF as in Theory, they should be the same, its the only explenation I have
#for why the splineFOAM isnt working.
import numpy as np
import imageTools as IT

splineDir='/home/stephen/OpenFOAM/Simulations2/Rumens/SplineChecks/RumenSplineTestCheck3/constant/patchDisplacements/'
positionDir='/home/stephen/OpenFOAM/Simulations2/Rumens/SplineChecks/RumenSplineTestCheck3/constant/patchDisplacements/'

patchName='inlet'

#read in the spline coeffs and calculte for node 0,0 at time 0,0.25,0.50,0.75 and 1
a=np.loadtxt(splineDir+patchName+'/a0')
b=np.loadtxt(splineDir+patchName+'/b0')
c=np.loadtxt(splineDir+patchName+'/c0')
d=np.loadtxt(splineDir+patchName+'/d0')

tRange=np.linspace(0,1,5)
xn=np.zeros(tRange.shape[0])

print a[0,0],b[0,0],c[0,0],d[0,0]
for i in range(tRange.shape[0]):
	xn[i]=a[0,0]*tRange[i]**3+b[0,0]*tRange[i]**2+c[0,0]*tRange[i]+d[0,0]

print tRange
print xn

#now we read in the two positions and do it that way
pos1=np.loadtxt(positionDir+patchName+'/patchDisplacement/0/patchDisplacement')
pos2=np.loadtxt(positionDir+patchName+'/patchDisplacement/1/patchDisplacement')
xPos=np.zeros(tRange.shape[0])
for i in range(tRange.shape[0]):
	xPos[i]=pos1[0,0]+((pos2[0,0]-pos1[0,0])*(tRange[i]))

print xPos
print pos2[0,0]-pos1[0,0]

print '++++++++++++++Error Marix++++++++++++'
t=0.5
errorMatrix=np.ones((a.shape[0],a.shape[1]))
for i in range(a.shape[0]):
	for j in range(a.shape[1]):
		sV=a[i,j]*t**3+b[i,j]*t**2+c[i,j]*t+d[i,j]
		pV=pos1[i,j]+((pos2[i,j]-pos1[i,j])*t)
		#print sV,pV
		errorMatrix[i,j]=(sV-pV)

	#print errorMatrix[i,:]

print np.sum(errorMatrix)

