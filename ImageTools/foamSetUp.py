#script to just convert stuff on the fly
import FoamTools as FT
workingDir='/home/stephen/OpenFOAM/Simulations2/Rumens/RumenSplineTestCheck/constant/patchDisplacements'
FT.writeXYZtoPointsVectorField(workingDir,'inletDisplacement','inletPositionNew')
FT.writeXYZtoPointsVectorField(workingDir,'wallDisplacement','wallPositionNew')