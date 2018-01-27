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
region='CranialSac'

#The number of sub points between significant motility bits, like max contraction, max expansion etc
numP=5
nSections=4 #The number of contraction events that occur

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
E1Times=[]
E2Times=[]
C1Times=[]
EndTimes=[]
numContours=[]
trajXList=[]
trajYList=[]

if SanityPlots:
	f,ax1=plt.subplots(1,1)

if SanityAnimation:
	f2,ax2=plt.subplots(1,1)

#loop through each scan and normalise everything to the same scale
for line in activeContours:
	print 'Opening contour:',line.rstrip()
	contractionList=open(workingDir+region+'/'+line.rstrip()+'/ContractionList')
	#loop through each contraction within each scan
	for contraction in contractionList: 
		print 'opening',contraction.rstrip()
		fid=open(workingDir+region+'/'+line.rstrip()+'/'+line.rstrip()+'_'+contraction.rstrip())
		count=0
		#get the zPos and Contraction length
		Pos=fid.readline().split()[2]
		numContours.append(fid.readline().split()[2])

		#set up vectors for storing the area
		areas=np.zeros(int(numContours[-1]))

		for fidLine in fid:
			if fidLine[0]=="#":
				temp=fidLine.split()
				#print np.shape(temp)
				areas[count]=float(temp[8])
				count=count+1

		fid.close()

		#First we find the min
		C1Time=np.argmin(areas[1:])+1
		
		#Then we find the max of each side of the min point
		E1Time=np.argmax(areas[:C1Time])
		E2Time=np.argmax(areas[C1Time:])+C1Time
		EndTime=count-1

		#print E1Time,C1Time,E2Time	
		C1Times.append(C1Time)
		E1Times.append(E1Time)
		E2Times.append(E2Time)
		EndTimes.append(EndTime)

		if (E2Time/(count-1))>=1:
			print 'Contraction',line,'second expansion ends and time end'

		zPos,dump,contours=IT.readContourFile(workingDir+region+'/'+line.rstrip()+'/'+line.rstrip()+'_'+contraction.rstrip())
		#Rest.append(contours[0])
		#C1Times.append(contours[C1Time])
		#E1Times.append(contours[E1Time])
		#E2Times.append(contours[E2Time])

		if SanityPlots:
			#ax1.plot(np.asarray(range(count))/float(count-1),(areas/areas[0]),'-')
			#ax1.plot(C1Time/float(count-1),areas[C1Time]/areas[0],'or')
			#ax1.plot(E1Time/float(count-1),areas[E1Time]/areas[0],'ob')
			#ax1.plot(E2Time/float(count-1),areas[E2Time]/areas[0],'ok')
			ax1.plot(np.asarray(range(count)),(areas/areas[0]),'-')
			ax1.plot(C1Time,areas[C1Time]/areas[0],'or')
			ax1.plot(E1Time,areas[E1Time]/areas[0],'ob')
			ax1.plot(E2Time,areas[E2Time]/areas[0],'ok')

		#Now we normalise data for each section of the contraction.

		#Tricky part. interp the the path of each time point for each normalised section.

		#start by getting the number of points in a contour
		cShape=np.shape(contours)
		#set the number of seed points, this is pretty arbitrary atm
		print 'Contour dims:'
		print cShape

		#interp each point in each contour accross time
		trajX=np.zeros((cShape[1],(nSections*numP)+1))
		trajY=np.zeros((cShape[1],(nSections*numP)+1))

		for node in range(cShape[1]):
			xMapFunc=interp1d(range(len(contours)),contours[:,node,0],kind='cubic')
			yMapFunc=interp1d(range(len(contours)),contours[:,node,1],kind='cubic')

			#Now we subdevide each section, could loop it, but since its location specific
			#I will hardcode this at least the first time, might fix up if I have time
			#later
			#rest positions
			trajX[node,0]=contours[0,node,0]
			trajY[node,0]=contours[0,node,1]

			#rest to expansion 1
			tPoints=np.linspace(0,E1Time,numP+1)
			xE1Points=xMapFunc(tPoints[1:])
			yE1Points=yMapFunc(tPoints[1:])

			#expansion 1->contraction 1
			tPoints=np.linspace(E1Time,C1Time,numP+1)
			xC1Points=xMapFunc(tPoints[1:])
			yC1Points=yMapFunc(tPoints[1:])

			#contraction 1-> expansion 2
			tPoints=np.linspace(C1Time,E2Time,numP+1)
			xE2Points=xMapFunc(tPoints[1:])
			yE2Points=yMapFunc(tPoints[1:])

			#expansion 2 -> end
			tPoints=np.linspace(E2Time,count-1,numP+1)
			xEndPoints=xMapFunc(tPoints[1:])
			yEndPoints=yMapFunc(tPoints[1:])

			#now put it all into the trajectory array, I know this could repace the Points stuff above....
			trajX[node,1:numP+1]=xE1Points
			trajX[node,numP+1:(numP*2)+1]=xC1Points
			trajX[node,(numP*2)+1:(numP*3)+1]=xE2Points
			trajX[node,(numP*3)+1:(numP*4)+1]=xEndPoints
			
			trajY[node,1:numP+1]=yE1Points
			trajY[node,numP+1:(numP*2)+1]=yC1Points
			trajY[node,(numP*2)+1:(numP*3)+1]=yE2Points
			trajY[node,(numP*3)+1:(numP*4)+1]=yEndPoints

		#Now smooth it all into a list for each contraction
		trajXList.append(trajX)
		trajYList.append(trajY)
		if SanityAnimation:
			plt.ion()
			f2.suptitle('Contraction '+line+' '+contraction)
			for t in range(numP*nSections):
				ax2.clear()
				ax2.set_xlim(50,350)
				ax2.set_ylim(-400,-100)
				ax2.plot(trajX[:,0],trajY[:,0]*-1,'rx')
				ax2.plot(trajX[:,t],trajY[:,t]*-1,'ob')
				plt.pause(0.2)

			plt.pause(2)
   			plt.ioff()

if SanityAnimation:
	#lets see it all smooshed togeather in glory
	plt.ion()
	f2.suptitle('IT LIVES!!!!')
	for t in range(numP*nSections):
		ax2.clear()
		ax2.set_xlim(50,350)
		ax2.set_ylim(-400,-100)
		for i in range(np.shape(trajXList)[0]):
			ax2.plot(trajXList[i][:,0],trajYList[i][:,0]*-1,'rx')
			ax2.plot(trajXList[i][:,t],trajYList[i][:,t]*-1,'ob')
		plt.pause(0.2)

#Lets have some stats
if SampleStats:
	print '+++Sample stats for Cranial Sac+++\n',
	print 'average time to first expansion:',np.mean(E1Times), np.std(E1Times)
	print 'average time to first contraction:',np.mean(C1Times), np.std(C1Times)
	print 'average time to second expansion:',np.mean(E2Times), np.std(E2Times)
	print 'average time to end of sequence:',np.mean(EndTimes), np.std(EndTimes)


#This section is used to some sanity PSLR stuff.
SanityPSLR=True
location=(numP)-2 #contraction peak

if SanityPSLR:

	#Constract a set of rest shapes, and contraction shapes
	cRest=[0]*np.shape(trajXList)[0]
	cCont=[0]*np.shape(trajXList)[0]
	f3,ax3=plt.subplots(2,2)
	for i in range(np.shape(trajXList)[0]):
		temp=np.zeros((np.shape(trajXList)[1],2))
		temp[:,0]=trajXList[i][:,0]
		temp[:,1]=trajYList[i][:,0]
		cRest[i]=temp

		temp=np.zeros((np.shape(trajXList)[1],2))
		temp[:,0]=trajXList[i][:,location]
		temp[:,1]=trajYList[i][:,location]
		cCont[i]=temp
		#plot the contracted state
		ax3[0,0].plot(cRest[i][:,0],cRest[i][:,1]*-1,'rx')
		ax3[0,0].plot(cCont[i][:,0],cCont[i][:,1]*-1,'ok')
		#ax3_1.set_xlim(100,300)
		#ax3_1.set_ylim(-350,-150)
		ax3[0,0].axis('equal')

	#Now align all the rest starts
	cShape=np.shape(cRest)
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

	for i in range(cShape[0]):
		ax3[0,1].plot(rPoints[i][:,0],rPoints[i][:,1]*-1,'-or')
		ax3[0,1].plot(cPoints[i][:,0],cPoints[i][:,1]*-1,'-xb')
		ax3[0,1].axis('equal')

	#Perform PSLR on aligned geometries

	#Step 1: Flatten the input vectors
	t0=[0]*(cShape[0]-1)
	t1=[0]*(cShape[0]-1)

	TestShape=3
	count=0
	for i in range(cShape[0]):
		if i==TestShape:
			print 'skipping', i
		else:
			t0[count]=np.ravel(rPoints[i])
			t1[count]=np.ravel(cPoints[i])
			count+=1

	#Step 2: Perform PSL regression
	pls=PLSRegression(n_components=3)
	pls.fit(t0,t1)

	test=pls.predict(np.ravel(rPoints[TestShape]).reshape(1,-1)) #PSLR needs to work with vectors in teh form [[]], its weird
	test=test.reshape((cShape[1],cShape[2]))

	ax3[1,0].plot(rPoints[TestShape][:,0],rPoints[TestShape][:,1]*-1,'ob')
	ax3[1,0].plot(cPoints[TestShape][:,0],cPoints[TestShape][:,1]*-1,'or')
	ax3[1,0].plot(test[:,0],test[:,1]*-1,'xk')
	ax3[1,0].axis('equal')

	#Calculate error 
	RMSE=np.sqrt(np.mean((cPoints[TestShape]-test)**2))
	CS=IT.calculateCentroidSize(cPoints[TestShape][:,0],cPoints[TestShape][:,1])
	print RMSE
	#print CS
	#print RMSE/CS
	#Calcualte the bounding box
	bbox=[np.min(cPoints[TestShape][:,0]),np.max(cPoints[TestShape][:,0]),np.min(cPoints[TestShape][:,1]),np.max(cPoints[TestShape][:,1])]
	unitLength=np.sqrt((bbox[0]-bbox[1])**2 + (bbox[2]-bbox[3])**2)
	print unitLength
	print RMSE/unitLength*100,'%'




#Now we apply PSLR to each time segment


plt.show()




		 






