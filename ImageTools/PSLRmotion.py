#This script is to try and see if we can interpolate motion from cine-CTs using Partial Least Squares regression.
import numpy as np
from sklearn.cross_decomposition import PLSRegression
import imageTools as IT
import matplotlib.pyplot as plt

#Read in the reticulum files, detect the zero position and maximum contraction.

workingDir = '/home/stephen/Documents/PhD/Data/Motility/SortedByLocation/'
region='Reticulum'

activeContours=open(workingDir+region+'/ProcessList').readlines()

peakTime=[]
areasAll=[]
cRest=[]
cCont=[]

for line in activeContours:
	print 'Opening contour:',line.rstrip()
	contractionList=open(workingDir+region+'/'+line.rstrip()+'/ContractionList')
	for contraction in contractionList: 
		print 'opening',contraction.rstrip()
		fid=open(workingDir+region+'/'+line.rstrip()+'/'+line.rstrip()+'_'+contraction.rstrip())
		count=0
		#get the zPos and Contraction length
		Pos=fid.readline().split()[2]
		numContours=fid.readline().split()[2]
		areas=np.zeros(int(numContours))
		CS=np.zeros(int(numContours))
		for fidLine in fid:
			if fidLine[0]=="#":
				temp=fidLine.split()
				#print np.shape(temp)
				areas[count]=float(temp[8])
				CS[count]=float(temp[6])
				count=count+1

		fid.close()

		peakTime.append(np.argmin(areas))
		areasAll.append(areas)

		zPos,numContours,contours=IT.readContourFile(workingDir+region+'/'+line.rstrip()+'/'+line.rstrip()+'_'+contraction.rstrip())
		cRest.append(contours[0])
		cCont.append(contours[np.argmin(areas)])

print np.size(cRest),np.shape(cCont)

f,ax1=plt.subplots(1,1)
flib=8
ax1.plot(cRest[flib][:,0],cRest[flib][:,1]*-1,'ob')
ax1.plot(cCont[flib][:,0],cCont[flib][:,1]*-1,'or')

#Step one, flatten cRest to 1D vector x1,y1,z1,x2,y2,z2 and same for cCont
cShape=np.shape(cRest)
#t0=np.reshape(cRest[0],(1,cShape[1]*cShape[2]))
t0=[0]*(cShape[0])
t1=[0]*(cShape[0])

for i in range(cShape[0]):
	#t0[i]=np.reshape(cRest[i],(1,cShape[1]*cShape[2]))
	t0[i]=np.ravel(cRest[i])
	t1[i]=np.ravel(cCont[i])

print np.shape(t0)

pls2 = PLSRegression(n_components=2)
pls2.fit(t0, t1)

test=pls2.predict(np.ravel(cRest[flib]).reshape(1,-1))
testPlot=test.reshape((cShape[1],cShape[2]))
ax1.plot(testPlot[:,0],testPlot[:,1]*-1,'xk')

#plot some like vector lines shit 
#for i in range(cShape[1]):
	#ax1.plot([testPlot[i,0],cCont[flib][i,0]],[testPlot[i,1]*-1,cCont[flib][i,1]*-1],'--k')

#print t0.reshape((cShape[1],cShape[2]))



plt.show()
