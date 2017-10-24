#Scrip is used to process motility files to look at changing area and what not.
import matplotlib.pyplot as plt
import numpy as np
import imageTools as IT

workingDir = '/home/stephen/Documents/PhD/Data/Motility/SortedByLocation/'
#pick your region
region='Reticulum'
#open the Processlist

activeContours=open(workingDir+region+'/ProcessList').readlines()
legendString=[]

for line in activeContours:
	print 'Opening contour:',line.rstrip()
	fid=open(workingDir+region+'/'+line.rstrip()+'/'+line.rstrip()+'_Contraction1')
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