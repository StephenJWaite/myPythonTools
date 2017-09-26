#script writes a xyz file to a FOAM point vector field
import numpy as np
print 'running writeXYZtoPointVectorField...'
workingDir='/home/stephen/OpenFOAM/Simulations2/Rumens/NewRumenNesiSetUp/constant/patchDisplacements'
filename='wall'
#read in the xyz file
xyzPoints=np.loadtxt(workingDir+'/'+filename+'Displacement')

#OpenFOAM file
fid=open(workingDir+'/'+filename+'PositionNew','w')
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




