#script replaces xyz coords in ply file
import numpy as np
workingDir='/home/stephen/Documents/Python/myPythonTools/Geoms/SoblovTests/'
filenamePLY='FixAttempt2'
filenameXYZ='BMES3'
fid=open(workingDir+filenamePLY+'.ply','r')
fidwrite=open(workingDir+filenameXYZ+'_outfile.ply','w')
flag=0
count=0
#loop through line by line
for line in fid:
	if flag==1:
		linetemp=map(float,line.rstrip().split())
		#print np.shape(temp)
		if np.shape(linetemp)[0]==3:
			temp=vertexArray[count,:]
			fidwrite.write('{0} '.format(temp[0]))
			fidwrite.write('{0} '.format(temp[1]))
			fidwrite.write('{0}\n'.format(temp[2]))
			count=count+1
		else:
			flag=2
			fidwrite.write(line)
	else:
		fidwrite.write(line)
			
		#we move values into an array
	if line[0:14]=='element vertex':
		vertexNum=int(line[14:].rstrip())
		#vertexArray=np.zeros((vertexNum,3))
		#read in xyz file
		vertexArray=np.loadtxt(workingDir+filenameXYZ+'.xyz')
		#sanit check dimensions
		print 'vertexNum:',vertexNum,'vertexArray:',vertexArray.shape

	if line.rstrip()=='end_header':
		print 'hello'
		flag=1

	#here we write to the new file
	

#export xyz
print 'finished reading'
