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

#Step 1: Read in the experimental data
workingDir = '/home/stephen/Documents/PhD/Data/Motility/SortedByLocation/'
region='CranialSac'

SanityPlots=True

activeContours=open(workingDir+region+'/ProcessList').readlines()
print 'Running fitting cranial sac motion with PSLR for:'
for line in activeContours:
	print line.rstrip()

Rest=[]
E1Times=[]
E2Times=[]
C1Times=[]

if SanityPlots:
	f,ax1=plt.subplots(1,1)

#loop through each scan
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
		numContours=fid.readline().split()[2]

		#set up vectors for storing the area
		areas=np.zeros(int(numContours))

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

		print E1Time,C1Time,E2Time

		
		C1Times.append(C1Time)
		E1Times.append(E1Time)
		E2Times.append(E2Time)

		zPos,numContours,contours=IT.readContourFile(workingDir+region+'/'+line.rstrip()+'/'+line.rstrip()+'_'+contraction.rstrip())
		Rest.append(contours[0])
		C1Times.append(contours[C1Time])
		E1Times.append(contours[E1Time])
		E2Times.append(contours[E2Time])

		if SanityPlots:
			ax1.plot(np.asarray(range(count))/float(count-1),(areas/areas[0]),'-')
			ax1.plot(C1Time/float(count-1),areas[C1Time]/areas[0],'or')
			ax1.plot(E1Time/float(count-1),areas[E1Time]/areas[0],'ob')
			ax1.plot(E2Time/float(count-1),areas[E2Time]/areas[0],'ok')

plt.show()




