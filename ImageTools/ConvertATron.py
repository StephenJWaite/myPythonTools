#ConvertATron is a script to convert some xyzs to vectorFields.
import FoamTools as FT
#workingDir='/home/stephen/OpenFOAM/Simulations2/Rumens/RumenSplineTest/constant/manualSplineTests/'
#workingDir='/home/stephen/OpenFOAM/Simulations2/Rumens/RumenSplineTest/constant/patchPositions/'
#workingDir='/home/stephen/OpenFOAM/Simulations2/Rumens/SplineChecks/RumenSplineTestCheck3/constant/patchDisplacements/'
workingDir='/home/stephen/OpenFOAM/Simulations2/Rumens/SplineChecks/RumenSplineTestCheck6/constant/patchDisplacements/'

#FT.writeXYZtoPointsVectorField(workingDir+'inlet/patchDisplacement/1','patchDisplacement','inletPositionNew')
#FT.writeXYZtoPointsVectorField(workingDir+'wall/1','patchDisplacements','wallPositionNew')

FT.writeXYZtoPointsVectorField(workingDir+'wall/0','a','wallA')
FT.writeXYZtoPointsVectorField(workingDir+'wall/0','b','wallB')
FT.writeXYZtoPointsVectorField(workingDir+'wall/0','c','wallC')
FT.writeXYZtoPointsVectorField(workingDir+'wall/0','d','wallD')

FT.writeXYZtoPointsVectorField(workingDir+'inlet/0','a','inletA')
FT.writeXYZtoPointsVectorField(workingDir+'inlet/0','b','inletB')
FT.writeXYZtoPointsVectorField(workingDir+'inlet/0','c','inletC')
FT.writeXYZtoPointsVectorField(workingDir+'inlet/0','d','inletD')

#FT.writeXYZtoPointsVectorField(workingDir,'inletDisplacement','inletPositionNew')
#FT.writeXYZtoPointsVectorField(workingDir,'wallDisplacement','wallPositionNew')