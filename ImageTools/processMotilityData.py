#Scrip is used to process motility files to look at changing area and what not.
import matplotlib.pyplot as plt
import numpy as np
import imageTools as IT

workingDir = '/home/stephen/Documents/PhD/Data/Motility/SortedByLocation/'
#pick your region
region='CranialSac'
#open the Processlist

f,ax1=plt.subplots(1,1)

activeContours=open(workingDir+region+'/ProcessList').readlines()
print activeContours
legendString=[]
meanReduction=[]
reductionLocation=[]
peakTime=[]
relaxTime=[]
areasAll=[]
CSAll=[]
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

		ax1.plot(np.asarray(range(count))/float(count-1),(areas/areas[0]),'-')
		legendString.append(line.rstrip())
		fid.close()

		#Calculate the mean area reduction
		maxValue=np.max(areas)/areas[0]
		minValue=np.min(areas)/areas[0]

		if (np.abs(maxValue-1) > np.abs(1-minValue)):
			meanReduction.append(maxValue)
			reductionLocation.append(np.argmax(areas)/float(count-1))
			#ax1.plot(temp/float(count-1),maxValue,'ro')
		else:
			#temp=np.argmin(areas)
			reductionLocation.append(np.argmin(areas)/float(count-1))
			meanReduction.append(minValue)
			#ax1.plot(temp/float(count-1),minValue,'bo')


		print line.rstrip(), meanReduction[-1]
		peakTime.append(np.argmin(areas))
		relaxTime.append(count-peakTime[-1])
		areasAll.append(areas)
		CSAll.append(CS)


ax1.legend(legendString)

#ax1.plot(reductionLocation,meanReduction,'ob')

print 'average magnitude:', np.mean(meanReduction), np.std(meanReduction)
print 'average peak contraction time:', np.mean(peakTime), np.std(peakTime)
print 'average relaxation time:', np.mean(relaxTime), np.std(relaxTime)

f2,ax3=plt.subplots(1,1)
flib=4
ax3.plot(range(len(areasAll[flib])),areasAll[flib]/areasAll[flib][0],'-k')
ax3.plot(range(len(areasAll[flib])),CSAll[flib]/CSAll[flib][0],'-r')

ax3.legend(['Area','Centroid Size'])
plt.ylabel('area (mm2')
plt.xlabel('time (s)')
plt.title(legendString[flib])



print np.shape(areasAll)


plt.show()