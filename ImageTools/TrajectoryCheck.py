#Script is used to look at node trajectories and spline fits
import numpy as np
import imageTools as IT
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy import poly1d
from scipy.interpolate import interp1d

#Set working directory.
workingDir='/home/stephen/OpenFOAM/Simulations2/Rumens/RumenBMES/constant/patchDisplacements'
fileName='wall'
numDT=4
splineData=np.zeros((numDT))
xValue=0



for i in range(numDT):
	nodes=np.loadtxt(workingDir+'Tests/'+fileName+str(i))
	splineData[i]=nodes[xValue,0]
	plt.plot(i,nodes[xValue,0],'or')
	
#Calculate spline values
[a,b,c,d]=IT.basicCubicSpline(numDT,range(numDT),splineData)
IT.basicCubicSplinePlot(a,b.tolist(),c.tolist(),d.tolist(),range(numDT))


#fig3=plt.figure(3)
#ax3D=fig3.add_subplot(111, projection='3d')
#ax3D.plot(nodes[:,0],nodes[:,1],nodes[:,2],'.b')
#ax3D.plot(nodes[xValue:xValue+1,0],nodes[xValue:xValue+1,1],nodes[xValue:xValue+1,2],'or')
test=np.zeros(4)
#Now we are going to check the outfiles from me thingy
for i in range(numDT-2):
	#read in a b c and d from a time file
	a=np.loadtxt(workingDir+'CoeffsReduced/'+fileName+'/'+str(i)+'/a')
	b=np.loadtxt(workingDir+'CoeffsReduced/'+fileName+'/'+str(i)+'/b')
	c=np.loadtxt(workingDir+'CoeffsReduced/'+fileName+'/'+str(i)+'/c')
	d=np.loadtxt(workingDir+'CoeffsReduced/'+fileName+'/'+str(i)+'/d')

	aC=np.loadtxt(workingDir+'TransformedReduced/'+fileName+'/'+str(i)+'/a')
	bC=np.loadtxt(workingDir+'TransformedReduced/'+fileName+'/'+str(i)+'/b')
	cC=np.loadtxt(workingDir+'TransformedReduced/'+fileName+'/'+str(i)+'/c')
	dC=np.loadtxt(workingDir+'TransformedReduced/'+fileName+'/'+str(i)+'/d')

	#now MATH!
	numTime=20
	dt=np.linspace(i,i+1,numTime)
	xTraj=np.zeros(numTime)
	for k in range(numTime):
		xTraj[k]=d[xValue,0]*(dt[k]**3)+c[xValue,0]*(dt[k]**2)+d[xValue,0]*dt[k]+a[xValue,0]

	#other check
	print 'i:',i
	root = poly1d(i,True)
	print root
	poly = 0
	print poly
	poly = d[xValue,0].tolist()*(root)**3
	print poly
	poly = poly + c[xValue,0].tolist()*(root)**2
	print poly
	poly = poly + b[xValue,0].tolist()*root
	print poly
	poly = poly + a[xValue,0]
	print poly
	print poly.c
	test[0],test[1],test[2],test[3]=poly.c
	print test
	print 'aC:',aC[xValue,0],'bC:',bC[xValue,0],'cC:',cC[xValue,0],'dC:',dC[xValue,0]
	[aCC,bCC,cCC,dCC]=IT.transformCubicCoeffiencets(a[xValue,:2],b[xValue,:2],c[xValue,:3],d[xValue,:2],dt[:2])
	print 'aCC:',aCC,'bCC:',bCC,'cCC:',cCC,'dCC:',dCC
	#print 'poly:',poly
	#set up a vector of x values from xi to xi+1
	x_step=0.05
	x_pts = np.arange(i,(i+1)+x_step,x_step)
	plt.plot(x_pts,poly(x_pts),'xb')
	#print 'a:',a[xValue,0],'b:',b[xValue,0],'c:',c[xValue,0],'d:',d[xValue,0]

	print poly

	#plt.plot(dt,xTraj,'xb')

plt.show()