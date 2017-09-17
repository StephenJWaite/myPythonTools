#Script for extracting the xyz coords from a ply file.
#Read in a ply file
import numpy as np
workingDir='/home/stephen/Documents/Python/myPythonTools/imageToolsTestSpace/Geoms/FittingResults/'
filename='1_2_211'
fid=open(workingDir+filename+'.ply')
flag=0
count=0
#loop through line by line
for line in fid:
	if flag==1:
		temp=map(float,line.rstrip().split())
		#print np.shape(temp)
		if np.shape(temp)[0]==3:
			vertexArray[count,:]=map(float,line.rstrip().split())
			count=count+1
		else:
			flag=2
			
		#we move values into an array
	if line[0:14]=='element vertex':
		vertexNum=int(line[14:].rstrip())
		vertexArray=np.zeros((vertexNum,3))
		print 'found ele vetrx',vertexNum, 'size:',vertexArray.shape
	if line.rstrip()=='end_header':
		print 'hello'
		flag=1

#export xyz
np.savetxt(workingDir+filename+'_outfile.xyz',vertexArray)
print 'finished reading'
