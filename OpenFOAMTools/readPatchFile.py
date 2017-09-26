#Function is used to read a specified patch file, and export the xyz coords Read patch file
import numpy as np
print 'running readPatchFile'
workingDir='/home/stephen/OpenFOAM/Simulations2/Rumens/NewRumen/constant'
filename='wallPatchNodes'
fid=open(workingDir+'/'+filename)
flag=0
count=0

for line in fid:

	if flag==1 and line.rstrip()==")":
			flag=2

	if flag==1:
		if line.strip("(").rstrip():
			#print line.rstrip().strip("()").split()
			#print np.size(map(float,line.rstrip().strip("()").split()))
			patchArray[count,:]=map(float,line.rstrip().strip("()").split())
			#rint patchArray[count,:]
			count=count+1

		

	if line.rstrip().isdigit():
		arraySize=int(line)
		print 'Patch has', arraySize, 'nodes'
		patchArray=np.zeros((arraySize,3))
		flag=1
	
#print out xyz
np.savetxt(workingDir+'/'+filename+'_outfile.xyz',patchArray)
print 'finished reading'
#print 'Patch array is ',np.shape(patchArray)
#print patchArray