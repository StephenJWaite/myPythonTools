#script reads in centroid sizes and areas over time, and plots them.
import matplotlib.pyplot as plt
import numpy as np
import imageTools as IT
#set workingDir
workingDir = './../imageToolsTestSpace/'

activeContours=open(workingDir+'ContourList').readlines()
legendString=[]

f,(ax1,ax2)=plt.subplots(2,1)

for line in activeContours:
	print 'Opening contour:',line.rstrip()
	fid=open(workingDir+line.rstrip()+'/Results.txt')
	count=0
	#get the zPos and Contraction length
	Pos=fid.readline().split()[2]
	numContours=fid.readline().split()[2]
	areas=np.zeros(int(numContours))
	for fidLine in fid:
		if fidLine[0]=="#":
			temp=fidLine.split()
			#print np.shape(temp)
			areas[count]=float(temp[8])
			count=count+1
	print count
	ax1.plot(range(count),areas,'-')
	print 'moooooooo'
	print line
	print legendString
	legendString.append(line.rstrip())
	print legendString
	fid.close()
	#zPos,numContours,contours=IT.readContourFile(workingDir+line.rstrip()+'/Results.txt')
	#store the zPos and numContous in contourInformation
	#contourInformation=contourInformation+[[zPos,numContours]]
	#contourData=contourData+[contours]

ax1.legend(legendString)

#now we will plot how that looks with the contraction map
CMap=np.loadtxt(workingDir+'contractionMap')
#set up infromation arrays
contourData=[]
contourInformation=[]

for line in activeContours:
	print 'Opening contour:',line.rstrip()
	zPos,numContours,contours=IT.readContourFile(workingDir+line.rstrip()+'/Results.txt')
	#store the zPos and numContous in contourInformation
	contourInformation=contourInformation+[[zPos,numContours]]
	contourData=contourData+[contours]

contractionMatrix,endtime=IT.createContractionMatrix(activeContours,contourInformation,contourData,CMap)
Cshape=np.shape(contractionMatrix)

#and now we calculate the vector of areas for each contour, and plot.
for i in range(Cshape[1]):
	areas=np.zeros(Cshape[0])
	for j in range(Cshape[0]):
		areas[j]=IT.calculateClosedContourArea(contractionMatrix[j][i][:,0],contractionMatrix[j][i][:,1])

	ax2.plot(range(Cshape[0]),areas,'-')


plt.show()