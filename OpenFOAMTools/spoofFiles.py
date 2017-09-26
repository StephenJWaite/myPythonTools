#File is used to test static wall motion....it spoofs files
import numpy as np
print 'running writeXYZtoPointVectorField...'
workingDir='/home/stephen/OpenFOAM/Simulations2/Rumens/NewRumen/constant/patchDisplacements'
filename='wall'
#read in the xyz file
xyzPoints=np.loadtxt(workingDir+'/'+filename+'Displacement')

#spoofing a OpenFOAM file
fid=open(workingDir+'/'+filename+'spoofFoam','w')
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
	temp='0.0 0.0 0.0'	
	fid.write('('+temp+')\n')

fid.write(')')

#spoof a point dispalcement file
pDisp=np.zeros((xyzPoints.shape[0],3))
np.savetxt(workingDir+'/'+filename+'spoofpDisp',pDisp)