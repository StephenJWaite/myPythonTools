#Script uses paraview simple to convert a set of ply files numbered 0->n to VTK files for animatiing wall motion
#Note, you must invoke a OF alias first of40 to use pvpython
from paraview.simple import *
import numpy as np

workingDir='/home/stephen/OpenFOAM/Simulations2/MotionCheck/OutFiles/'

numFiles=15

for i in range(0,numFiles):
	print 'Processing image:',i
	plyFile=PLYReader(FileName=str(i)+'.ply')
	writer=XMLPolyDataWriter(plyFile,FileName='test_'+str(i)+'.vtp')
	writer.UpdatePipeline()