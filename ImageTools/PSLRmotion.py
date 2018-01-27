#This script is to try and see if we can interpolate motion from cine-CTs using Partial Least Squares regression.
import numpy as np
from sklearn.cross_decomposition import PLSRegression
from sklearn import cross_validation
import imageTools as IT
import matplotlib.pyplot as plt

###############################################################################
#Read in the reticulum files, detect the zero position and maximum contraction.
workingDir = '/home/stephen/Documents/PhD/Data/Motility/SortedByLocation/'
region='Reticulum'

activeContours=open(workingDir+region+'/ProcessList').readlines()

peakTime=[]
areasAll=[]
cRest=[]
cCont=[]
legend=[]

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

		legend.append(line.rstrip())
		peakTime.append(np.argmin(areas))
		areasAll.append(areas)

		zPos,numContours,contours=IT.readContourFile(workingDir+region+'/'+line.rstrip()+'/'+line.rstrip()+'_'+contraction.rstrip())
		cRest.append(contours[0])
		cCont.append(contours[np.argmin(areas)])

print np.size(cRest),np.shape(cCont)

f,ax1=plt.subplots(1,1)
flib=1

#Going to try flipping 6905_07 -> [1:3]
#it flipped it back, cause of the landmark numbering lol



ax1.plot(cRest[flib][:,0],cRest[flib][:,1]*-1,'ob')
ax1.plot(cCont[flib][:,0],cCont[flib][:,1]*-1,'or')
#############################################################################

#############################################################################
#Alignment procedures of rest and contraction position





#Step 1: set up storeage vectors for rest and contraction information, and 
#        perform 2D alignment.

cShape=np.shape(cRest)

#Rest state
rR=[0]*cShape[0]
rS=[0]*cShape[0]
rT=[0]*cShape[0]
rPoints=[0]*cShape[0]
rMeanShape=cRest[0] #start with the mean as the first shape

print 'aligning shapes at rest'
for i in range(cShape[0]):
	rPoints[i],rR[i],rS[i],rT[i]=IT.procrustesAlignment2D(rMeanShape,cRest[i])

#Contracted state
cR=[0]*cShape[0]
cS=[0]*cShape[0]
cT=[0]*cShape[0]
cPoints=[0]*cShape[0]
cMeanShape=cRest[0] 

print 'apply the alignment of rest to contractoin'
for i in range(cShape[0]):
	cPoints[i]=IT.apply2DTransformation(list(cCont[i]),rT[i],rR[i],rS[i])
#	cPoints[i],cR[i],cS[i],cT[i]=IT.procrustesAlignment2D(cMeanShape,cCont[i])


#plot the aligned shapes
f2,ax2=plt.subplots(1,1)

for i in range(cShape[0]):
	ax2.plot(rPoints[i][:,0],rPoints[i][:,1]*-1,'-or')
	ax2.plot(cPoints[i][:,0],cPoints[i][:,1]*-1,'-xb')

#############################################################################

#############################################################################
#Perform PSLR on aligned geometries

#Step 1: Flatten the input vectors
t0=[0]*(cShape[0]-1)
t1=[0]*(cShape[0]-1)

TestShape=2
count=0
for i in range(cShape[0]):
	if i==TestShape:
		print 'skipping', i,':',legend[i]
	else:
		t0[count]=np.ravel(rPoints[i])
		t1[count]=np.ravel(cPoints[i])
		count+=1

#Step 2: Perform PSL regression
pls=PLSRegression(n_components=2)
pls.fit(t0,t1)

#Step 3: Sanity check prediction of a random shape


samplePoints,sR,sS,sT=IT.procrustesAlignment2D(rMeanShape,cRest[TestShape])
test=pls.predict(np.ravel(samplePoints).reshape(1,-1)) #PSLR needs to work with vectors in the form [[]], its weird
test=test.reshape((cShape[1],cShape[2]))

#transform the real contraciton points
contractedPoints=IT.apply2DTransformation(cCont[TestShape],sT,sR,sS)

f3,ax3=plt.subplots(1,1)
ax3.plot(samplePoints[:,0],samplePoints[:,1]*-1,'ob')
ax3.plot(contractedPoints[:,0],contractedPoints[:,1]*-1,'or')
ax3.plot(test[:,0],test[:,1]*-1,'xk')
ax3.axis('equal')

#Calcualte the bounding box
bbox=[np.min(cPoints[TestShape][:,0]),np.max(cPoints[TestShape][:,0]),np.min(cPoints[TestShape][:,1]),np.max(cPoints[TestShape][:,1])]
RMSE=np.sqrt(np.mean((cPoints[TestShape]-test)**2))
unitLength=np.sqrt((bbox[0]-bbox[1])**2 + (bbox[2]-bbox[3])**2)
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


#n=cShape[0]
#kf_10 = cross_validation.KFold(n,n_folds=2, shuffle=False,random_state=1)
#mse=[]

#print 't0.shape:',np.shape(t0)
#print 't1.shape:',np.shape(t1)

#for i in np.arange(1,cShape[0]):
#	pls=PLSRegression(n_components=6)
#	scores=cross_validation.cross_val_score(pls,t0,t1,cv=kf_10,scoring='r2')
#	print scores
#	mse.append(-scores)

#f4,ax4=plt.subplots(1,1)
#print np.array(mse).shape
#ax4.plot(np.arange(1,cShape[0]),np.array(mse),'-o')
#plt.xlabel('Number of components')
#plt.ylabel('MSE')
#plt.title('asfsdf')







plt.show()
