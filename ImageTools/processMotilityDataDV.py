#Scrip is used to process motility files to look at changing area and what not.
import matplotlib.pyplot as plt
import numpy as np
import imageTools as IT

workingDir = '/home/stephen/Documents/PhD/Data/Motility/SortedByLocation/'
#pick your region
region='CentralRear'
#open the Processlist

f,(ax1,ax2)=plt.subplots(2,1)

activeContours=open(workingDir+region+'/ProcessListDV').readlines()
print activeContours
legendString=[]
meanReductionD=[]
peakTimeD=[]
relaxTimeD=[]
areasD=[]

meanReductionV=[]
peakTimeV=[]
relaxTimeV=[]
areasV=[]


for line in activeContours:
	print 'Opening contour:',line.rstrip()
	contractionList=open(workingDir+region+'/'+line.rstrip()+'/ContractionList')
	for contraction in contractionList: 
		#copy paste grossness
		print 'opening Dorsal',contraction.rstrip()
		fid=open(workingDir+region+'/'+line.rstrip()+'/'+line.rstrip()+'_'+contraction.rstrip()+'Dorsal')
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

		ax1.plot(np.asarray(range(count))/float(count-1),(areas/areas[0]),'-')
		legendString.append(line.rstrip())
		fid.close()

		#Calculate the mean area reduction
		meanReductionD.append(np.min(areas)/areas[0])
		peakTimeD.append(np.argmin(areas))
		relaxTimeD.append(count-peakTimeD[-1])
		areasD.append(areas)

		print 'opening Ventral',contraction.rstrip()
		fid=open(workingDir+region+'/'+line.rstrip()+'/'+line.rstrip()+'_'+contraction.rstrip()+'Ventral')
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

		ax2.plot(np.asarray(range(count))/float(count-1),(areas/areas[0]),'-')
		fid.close()

		#Calculate the mean area reduction
		meanReductionV.append(np.min(areas)/areas[0])
		peakTimeV.append(np.argmin(areas))
		relaxTimeV.append(count-peakTimeV[-1])
		areasV.append(areas)

ax1.legend(legendString)
ax2.legend(legendString)

print 'Dorsal Sats:'
print 'average magnitude:', np.mean(meanReductionD), np.std(meanReductionD)
print 'average peak contraction time:', np.mean(peakTimeD), np.std(peakTimeD)
print 'average relaxation time:', np.mean(relaxTimeD), np.std(relaxTimeD)

print 'Vental Sats:'
print 'average magnitude:', np.mean(meanReductionV), np.std(meanReductionV)
print 'average peak contraction time:', np.mean(peakTimeV), np.std(peakTimeV)
print 'average relaxation time:', np.mean(relaxTimeV), np.std(relaxTimeV)

f2,ax3=plt.subplots(1,1)
flib=4
ax3.plot(range(len(areasD[flib])),areasD[flib],'--k')
ax3.plot(range(len(areasV[flib])),areasV[flib],'-k')
ax3.legend(['Dorsal Sac','Ventral Sac'])
plt.ylabel('area (mm2')
plt.xlabel('time (s)')
plt.title(legendString[flib])


plt.show()