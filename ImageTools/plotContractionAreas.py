#script reads in centroid sizes and areas over time, and plots them.
import matplotlib.pyplot as plt
import numpy as np
import imageTools as IT
#set workingDir
workingDir = './../imageToolsTestSpace/'

activeContours=open(workingDir+'ContourList').readlines()

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
	plt.plot(range(count),areas,'-')
	fid.close()
	#zPos,numContours,contours=IT.readContourFile(workingDir+line.rstrip()+'/Results.txt')
	#store the zPos and numContous in contourInformation
	#contourInformation=contourInformation+[[zPos,numContours]]
	#contourData=contourData+[contours]


plt.show()