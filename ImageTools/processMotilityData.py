#Scrip is used to process motility files to look at changing area and what not.
import matplotlib.pyplot as plt
import numpy as np
import imageTools as IT

workingDir = '/home/stephen/Documents/PhD/Data/Motility/SortedByLocation/'
#pick your region
region='Reticulum'
#open the Processlist

f,ax1=plt.subplots(1,1)

activeContours=open(workingDir+region+'/ProcessList').readlines()
legendString=[]
meanReduction=[]
peakTime=[]
relaxTime=[]
contractionNum=0
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

	ax1.plot(range(count),(areas/areas[0]),'-')
	print line
	legendString.append(line.rstrip())
	fid.close()

	#Calculate the mean area reduction
	meanReduction.append(np.min(areas)/areas[0])
	peakTime.append(np.argmin(areas))
	relaxTime.append(count-peakTime[contractionNum])
	contractionNum=contractionNum+1


ax1.legend(legendString)

print 'average magnitude:', np.mean(meanReduction), np.std(meanReduction)
print 'average peak contraction time:', np.mean(peakTime), np.std(peakTime)
print 'average relaxation time:', np.mean(relaxTime), np.std(relaxTime)


plt.show()