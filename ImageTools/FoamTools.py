#FOAM tools are a series of functions to work with FOAM files in PYTHON (most file IO stuff)
import sys
import numpy as np

def readPatchFile(workingDir,filename):
	print 'running readPatchFile:'+filename
	fid=open(workingDir+filename)
	flag=0
	count=0
	patchArray=[]

	for line in fid:

		if flag==1 and line.rstrip()==")":
			flag=2

		if flag==1:
			if line.strip("(").rstrip():
				patchArray[count,:]=map(float,line.rstrip().strip("()").split())
				#increase count
				count=count+1

		if line.rstrip().isdigit():
			arraySize=int(line)
			print 'Patch has', arraySize, 'nodes'
			patchArray=np.zeros((arraySize,3))
			flag=1
	
	print 'finished reading'+filename
	return patchArray

def writeXYZtoPointsVectorField(workingDir,fileName,outFileName):
	#script writes a xyz file to a FOAM point vector field
	import numpy as np
	print 'running writeXYZtoPointVectorField...'

	#read in the xyz file
	xyzPoints=np.loadtxt(workingDir+'/'+fileName)

	#OpenFOAM file
	fid=open(workingDir+'/'+outFileName,'w')
	#Write header information
	fid.write('FoamFile\n')
	fid.write('{\n')
	fid.write('    version    2.0;\n')
	fid.write('    format     ascii;\n')
	fid.write('    class      vectorField;\n')
	fid.write('    location   "constant";\n')
	fid.write('    object     surfacePosition;\n')
	fid.write('}\n\n\n')
	#write the number of points
	fid.write('{0}\n'.format(xyzPoints.shape[0]))
	fid.write('(\n')
	for point in xyzPoints:
		temp=str(point).strip(' [] ')	
		fid.write('('+temp+')\n')

	fid.write(')')

def writeVectorField(workingDir,vectorField):
	print 'running writeVectorField...'
	#OpenFOAM file
	fid=open(workingDir,'w')
	#Write header information
	fid.write('FoamFile\n')
	fid.write('{\n')
	fid.write('    version    2.0;\n')
	fid.write('    format     ascii;\n')
	fid.write('    class      vectorField;\n')
	fid.write('    location   "constant";\n')
	fid.write('    object     surfacePosition;\n')
	fid.write('}\n\n\n')
	#write the number of points
	fid.write('{0}\n'.format(vectorField.shape[0]))
	fid.write('(\n')
	for point in vectorField:
		temp=str(point).strip(' [] ')	
		fid.write('('+temp+')\n')

	fid.write(')')

def spoofEmptyDisplacementFile(workingDir):
	print 'running spoofEmptyDisplacementFile <- someones get paralled ^^'
	fid=open(workingDir,'w')
	fid.write('{0}'.format(0))
	fid.close()


def createBlankVectorField(workingDir,filename):
	print 'running createBlankVectirField'

	#OpenFOAM file
	fid=open(workingDir+'/'+filename,'w')
	#Write header information
	fid.write('FoamFile\n')
	fid.write('{\n')
	fid.write('    version    2.0;\n')
	fid.write('    format     ascii;\n')
	fid.write('    class      vectorField;\n')
	fid.write('    location   "constant";\n')
	fid.write('    object     surfacePosition;\n')
	fid.write('}\n\n\n')
	#write the number of points
	fid.write('{0}()\n'.format(0))




	
