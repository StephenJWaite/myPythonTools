#Script for determing the motion of the cranial sac. The cranial sac goes through three motions
#Expansion
#Contraction
#Expansion
import sys
sys.path.insert(0,'./../../ImageTools')
import imageTools as IT
import numpy as np
from sklearn.cross_decomposition import PLSRegression
from sklearn import cross_validation
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import time

plt.switch_backend('TkAgg')

#Step 1: Read in the experimental data
workingDir = '/home/stephen/Documents/PhD/Data/Motility/SortedByLocation/'
region='CentralFront'

#The number of sub points between significant motility bits, like max contraction, max expansion etc
numP=5
nSections=2 #The number of contraction events that occur

#Turn this on if you want sanity plots
SanityPlots=True
SanityAnimation=False
SampleStats=True
SanityPSLR=True

activeContours=open(workingDir+region+'/ProcessList').readlines()
print 'Running fitting cranial sac motion with PSLR for:'
for line in activeContours:
	print line.rstrip()

Rest=[]
EndTimes=[]
numContours=[]
legendString=[]
cRest=[]
cCont=[]

if SanityPlots:
	f,ax1=plt.subplots(1,1)

if SanityAnimation:
	f2,ax2=plt.subplots(1,1)

#loop through each scan and normalise everything to the same scale
for line in activeContours:
	print 'Opening contour:',line.rstrip()
	contractionList=open(workingDir+region+'/'+line.rstrip()+'/ContractionList')
	#Open the file that says where each contractions max is
	cPointFid=open(workingDir+region+'/'+line.rstrip()+'/CPoint')

	#loop through each contraction within each scan
	for contraction in contractionList: 
		print 'opening',contraction.rstrip()
		zPos,dump,contours=IT.readContourFile(workingDir+region+'/'+line.rstrip()+'/'+line.rstrip()+'_'+contraction.rstrip())

		#Set the cPoint
		cPoint=int(cPointFid.readline().rstrip())

		EndTime=np.shape(contours)[0]

		EndTimes.append(EndTime)

		legendString.append(line.rstrip())

		#ax1.plot(contours[0][:,0],contours[0][:,1]*-1,'-o')
		cRest.append(contours[0])
		cCont.append(contours[cPoint])

		
#Now aligned the rest positions
cShape=np.shape(cRest)
rR=[0]*cShape[0]
rS=[0]*cShape[0]
rT=[0]*cShape[0]
rPoints=[0]*cShape[0]
cPoints=[0]*cShape[0]
rMeanShape=cRest[0] #start with the mean as the first shape
print 'aligning shapes at rest'
for i in range(cShape[0]):
	rPoints[i],rR[i],rS[i],rT[i]=IT.procrustesAlignment2D(rMeanShape,cRest[i])
	ax1.plot(rPoints[i][:,0],rPoints[i][:,1]*-1,'-ob')

	cPoints[i]=IT.apply2DTransformation(list(cCont[i]),rT[i],rR[i],rS[i])
	ax1.plot(cPoints[i][:,0],cPoints[i][:,1]*-1,'-or')

ax1.axis('equal')

#Step 1: Flatten the input vectors
t0=[0]*(cShape[0]-1)
t1=[0]*(cShape[0]-1)

TestShape=4
count=0
for i in range(cShape[0]):
	if i==TestShape:
		print 'skipping', i,':',legendString[i]
	else:
		t0[count]=np.ravel(rPoints[i])
		t1[count]=np.ravel(cPoints[i])
		count+=1

#Step 2: Perform PSL regression
pls=PLSRegression(n_components=5)
pls.fit(t0,t1)

#Step 3: Sanity check prediction of a random shape


samplePoints,sR,sS,sT=IT.procrustesAlignment2D(rMeanShape,cRest[TestShape])
test=pls.predict(np.ravel(samplePoints).reshape(1,-1)) #PSLR needs to work with vectors in the form [[]], its weird
test=test.reshape((cShape[1],cShape[2]))

#transform the real contraciton points
contractedPoints=IT.apply2DTransformation(cCont[TestShape],sT,sR,sS)

f3,ax3=plt.subplots(1,1)
ax3.plot(samplePoints[:,0],samplePoints[:,1]*-1,'-ob')
ax3.plot(contractedPoints[:,0],contractedPoints[:,1]*-1,'-or')
ax3.plot(test[:,0],test[:,1]*-1,'-xk')
ax3.axis('equal')

#Calcualte the bounding box
bbox=[np.min(cPoints[TestShape][:,0]),np.max(cPoints[TestShape][:,0]),np.min(cPoints[TestShape][:,1]),np.max(cPoints[TestShape][:,1])]
RMSE=np.sqrt(np.mean((cPoints[TestShape]-test)**2))
unitLength=np.sqrt((bbox[0]-bbox[1])**2 + (bbox[2]-bbox[3])**2)
print RMSE
print unitLength
print RMSE/unitLength*100,'%'

#Caluclate areas
areaR=IT.calculateClosedContourArea(rPoints[TestShape][:,0],rPoints[TestShape][:,1])
areaC=IT.calculateClosedContourArea(cPoints[TestShape][:,0],cPoints[TestShape][:,1])
areaT=IT.calculateClosedContourArea(test[:,0],test[:,1])

print 'AR', areaR
print 'AC:',areaC
print 'AT:',areaT
print 'Area percentrage:',(np.abs(areaC-areaT)/areaC)*100

		

plt.show()




		 






