#List of function tools for the sementor

import sys
import numpy as np
from numpy import poly1d #<- why is this different to numpy.poly1d?
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def helloWorld():
	print 'Hello Stephen, you nailed it'

def readPLY(meshString):
	#PLY have a certain format that we can take advantage of.
	#Step 1) Check its a ply by the header
	#Step 2) Read the number of nodes
	#Step 2) Read until the end hear line, we are now down with text
	f=open(meshString)

	#Check that its ply
	line=f.readline().rstrip()

	if line!='ply':
		print 'This is not a PLY file', line

		Nodes=0
		Elems=0
	else:
		print 'Reading mesh from path', meshString

		while line!='end_header':
			line=f.readline().rstrip()

			if line.find('element vertex')!=-1:
				nodeNum=int(line[15:])
				print 'Number of Nodes:', nodeNum

			if line.find('element face')!=-1:
				elemNum=int(line[13:])
				print 'Number of Elements:', elemNum

		print 'Reading Nodes...'
		Nodes=np.zeros(shape=(nodeNum,3))
		for x in range(0,nodeNum):
			line=f.readline().rstrip()
			temp=[float(s) for s in line.split()]
			Nodes[x,:]=temp

		print 'Reading Elements...'
		Elems=np.zeros(shape=(elemNum,3),dtype=np.int)
		for x in range(0,elemNum):
			line=f.readline().rstrip()
			temp=[int(s) for s in line.split()]
			Elems[x,:]=temp[1:]

		return Nodes, Elems

def extractPointsFromYPlane(Nodes,Elems,yPos):
	#This function is used to sort nodes on either side of
	#the yPos

	numNodes=Nodes.shape[0]
	numElems=Elems.shape[0]

	pointList=[]
	#loop through the elements
	for x in range(0,numElems):
		#These are surface tris
		temp=np.zeros(3)
		for y in range(0,3):
			if Nodes[Elems[x,y],1] < yPos:
				temp[y]= 1
		#Check to see if that temp array is all 1s (<Pos) or 0s (>Pos)
		#print temp, sum(temp)
		if (sum(temp) > 0) and (sum(temp) < 3):
			#print 'element spans the y plane'
			if sum(temp)==1:
				#we know that there is one point on the <yPos side
				#Note, argsort defaul goes smallest->biggest, so our
				#unique node is < Pos, ordered it will be the [2]
				#position in the elemOrdered
				elemOrdered=Elems[x,np.argsort(temp)]
				p1x, p1z = calculateYintercept(Nodes[elemOrdered[2],:],Nodes[elemOrdered[0],:],yPos)
				p2x, p2z = calculateYintercept(Nodes[elemOrdered[2],:],Nodes[elemOrdered[1],:],yPos)
				pointList.append([p1x,yPos,p1z])
				pointList.append([p2x,yPos,p2z])
			else:
				#Now our unique node is a 0, so its on the > Pos side.
				#When ordered it is now the [0] element in elemOrdered
				elemOrdered=Elems[x,np.argsort(temp)]
				p1x, p1z = calculateYintercept(Nodes[elemOrdered[0],:],Nodes[elemOrdered[2],:],yPos)
				p2x, p2z = calculateYintercept(Nodes[elemOrdered[0],:],Nodes[elemOrdered[1],:],yPos)
				pointList.append([p1x,yPos,p1z])
				pointList.append([p2x,yPos,p2z])
		else:
			pass#print 'No spanning'

	#Convert the list to a numpy array
	pointCloud=np.array(pointList)

	return pointCloud

def extractPointsFromZPlane(Nodes,Elems,zPos):
	#This function is used to sort nodes on either side of
	#the yPos

	numNodes=Nodes.shape[0]
	numElems=Elems.shape[0]

	pointList=[]
	#loop through the elements
	for x in range(0,numElems):
		#These are surface tris
		temp=np.zeros(3)
		for y in range(0,3):
			if Nodes[Elems[x,y],2] < zPos:
				temp[y]= 1
		#Check to see if that temp array is all 1s (<Pos) or 0s (>Pos)
		#print temp, sum(temp)
		if (sum(temp) > 0) and (sum(temp) < 3):
			#print 'element spans the y plane'
			if sum(temp)==1:
				#we know that there is one point on the <yPos side
				#Note, argsort defaul goes smallest->biggest, so our
				#unique node is < Pos, ordered it will be the [2]
				#position in the elemOrdered
				elemOrdered=Elems[x,np.argsort(temp)]
				#if isUnique(np.array([Nodes[elemOrdered[0],2],Nodes[elemOrdered[1],2],Nodes[elemOrdered[2],2]])):
				#	p1x, p1y = calculateZintercept(Nodes[elemOrdered[2],:],Nodes[elemOrdered[0],:],zPos,elemOrdered)
				#	p2x, p2y = calculateZintercept(Nodes[elemOrdered[2],:],Nodes[elemOrdered[1],:],zPos,elemOrdered)
				#	pointList.append([p1x,p1y,zPos])
				#	pointList.append([p2x,p2y,zPos])
				#else:
				#	print 'non Unique coordinate found removing...'
				p1x, p1y = calculateZintercept(Nodes[elemOrdered[2],:],Nodes[elemOrdered[0],:],zPos,elemOrdered)
				p2x, p2y = calculateZintercept(Nodes[elemOrdered[2],:],Nodes[elemOrdered[1],:],zPos,elemOrdered)
				pointList.append([p1x,p1y,zPos])
				pointList.append([p2x,p2y,zPos])

			else:
				#Now our unique node is a 0, so its on the > Pos side.
				#When ordered it is now the [0] element in elemOrdered
				elemOrdered=Elems[x,np.argsort(temp)]
				#if isUnique(np.array([Nodes[elemOrdered[0],2],Nodes[elemOrdered[1],2],Nodes[elemOrdered[2],2]])):
				#	p1x, p1y = calculateZintercept(Nodes[elemOrdered[0],:],Nodes[elemOrdered[2],:],zPos,elemOrdered)
				#	p2x, p2y = calculateZintercept(Nodes[elemOrdered[0],:],Nodes[elemOrdered[1],:],zPos,elemOrdered)
				#	pointList.append([p1x,p1y,zPos])
				#	pointList.append([p2x,p2y,zPos])
				#else:
				#	print 'non Unique coordinate found removing...'
				p1x, p1y = calculateZintercept(Nodes[elemOrdered[0],:],Nodes[elemOrdered[2],:],zPos,elemOrdered)
				p2x, p2y = calculateZintercept(Nodes[elemOrdered[0],:],Nodes[elemOrdered[1],:],zPos,elemOrdered)
				pointList.append([p1x,p1y,zPos])
				pointList.append([p2x,p2y,zPos])
		else:
			pass#print 'No spanning'

	#Convert the list to a numpy array
	pointCloud=np.array(pointList)

	return pointCloud


def isUnique(array):
	#Used to check if there are any double ups in the input array
	if array.size==np.unique(array).size:
		return True
	else:
		return False


def calculateYintercept(p1,p2,yPos):
	#Function is used to calculate the x and z coordinates of a y
	#intecept of two points

	#Start with xy plane
	xyTheta = np.arctan((p1[0]-p2[0])/(p1[1]-p2[1]))
	xPos=((yPos-p2[1])*np.tan(xyTheta))+p2[0]

	#zy plane
	zyTheta = np.arctan((p1[2]-p2[2])/(p1[1]-p2[1]))
	zPos=((yPos-p2[1])*np.tan(zyTheta))+p2[2]

	if np.isnan(xPos) or np.isnan(zPos):
		print 'Fucking NaN cat'
		print p1
		print p2
		print 'x values:', xPos, xyTheta, 'angle crap ', (p1[0]-p2[0])/(p1[1]-p2[1])
		print 'z values:', zPos, zyTheta

	return xPos, zPos

def calculateZintercept(p1,p2,zPos,elems):
	#Function is used to calculate the x and z coordinates of a y
	#intecept of two points

	#Start with xy plane
	xzTheta = np.arctan((p1[0]-p2[0])/(p1[2]-p2[2]))
	xPos=((zPos-p2[2])*np.tan(xzTheta))+p2[0]

	#zy plane
	yzTheta = np.arctan((p1[1]-p2[1])/(p1[2]-p2[2]))
	yPos=((zPos-p2[2])*np.tan(yzTheta))+p2[1]

	if np.isnan(xPos) or np.isnan(yPos):
		print 'Fucking NaN cat'
		print p1
		print p2
		print 'x values:', xPos, xzTheta, 'angle crap ', (p1[0]-p2[0])/(p1[1]-p2[1])
		print 'y values:', yPos, yzTheta


	if ((p1[2]-p2[2])==0.0):
	 	print 'DIV ZERO'
	 	print p2
		print 'x values:', xPos, xzTheta, 'angle crap ', (p1[0]-p2[0])/(p1[1]-p2[1])
		print 'y values:', yPos, yzTheta
		print 'elem numbers:', elems

	return xPos, yPos

#Function weightedGreyScale converts an RGB image into a weighted gray scale
def weightedGreyScale(imageData):
	#check that this is an RGB array.
	if imageData.shape[2]!=3:
		print 'Error in weightedGreyScale, image is not RGB'
		return -1
	#set up an empty 3D array
	grey=np.zeros((imageData.shape[0],imageData.shape[1]))
	#loop through the array x and y dimension and average
	for rownum in range(len(imageData)):
		for colnum in range(len(imageData[rownum])):
			grey[rownum][colnum] = 0.299*imageData[rownum][colnum] + 0.587*imageData[rownum][colnum] + 0.114*imageData[rownum][colnum]

	return grey

def pickPoint(fig):

	print 'running function pickPoint'

	def onclick(event):
		print '+Running onclick event+'
		print event.xdata, event.ydata
		fig.canvas.mpl_disconnect(cid)
		return event.xdata, event.ydata


	cid = fig.canvas.mpl_connect('button_press_event', onclick)
	return cid

#pickedPoints is a class for selecting and storing x,y coordinates from user selected locations on a figure.
#the class takes the number of points that the user wants to select and the figure handle for the figure
#they wish to select from.
#Still to implement, free form selection for unspecified number of points, breaks with a keystroke enter.
class pickedPoints():
	def __init__(self,numPoints,fig):
		#initialise self varaibles
		self.points=np.zeros((numPoints,2)) #array of x and y points
		self.count=0 #counter used to determine what point number the user is at.
		self.numPoints = numPoints# total number of points that the user will select
		self.fig=fig #Figure handle for the figure the user will select from
		self.freeSelect=False #no yet is use, true will be for free form selection
		#create button event for function pickPoints
		self.cid = self.fig.canvas.mpl_connect('button_press_event',self.pickPoints)

	def pickPoints(self,event):
		print 'Runing pickPoints from class pickedPoints'
		#assigned mouse click x and y values to point array
		self.points[self.count,:] = event.xdata,event.ydata
		#increase point counter
		self.count=self.count+1
		print event.xdata, event.ydata
		#excape condition, if point counter equals the number of points, disconnect button event
		if self.count==self.numPoints:
			self.fig.canvas.mpl_disconnect(self.cid)

class classTest():
	def __init__(self,test):
		self.cat=10/test
		print 'cat in class:',self.cat

def snapSphere(contourData,point,sR):
	#loop through the contour data, check if each point is inside the snap sphere
	print 'Running snapSphere...'
	sphereNodes=np.ones((1,2))*-1
	for i in range(len(contourData)):
		if (np.abs(contourData[i,0]-point[0])<=sR) and (np.abs(contourData[i,1]-point[1])<=sR):
			dist=np.sqrt((contourData[i,0]-point[0])**2 + (contourData[i,1]-point[1])**2)
			if sphereNodes[0,0]==-1:
				sphereNodes[0,:]=[dist,i]
			else:
				sphereNodes=np.append(sphereNodes,[[dist,i]],axis=0)

	#sort the sphereNodes
	sphereNodes=sphereNodes[np.argsort(sphereNodes[:,0]),:]
	#We take the closest points, and now only consider the nodes on either side. Special case if its at the end of the array
	if sphereNodes[0,1]==len(contourData)-1:
		print 'Contour Length',len(contourData)-1
		nodeGroup=np.asarray([sphereNodes[0,1],sphereNodes[0,1]-1,0]).astype(int)
	elif sphereNodes[0,1]==0:
		nodeGroup=np.asarray([sphereNodes[0,1],sphereNodes[0,1]+1,len(contourData)-1]).astype(int)
	else:
		nodeGroup=np.asarray([sphereNodes[0,1],sphereNodes[0,1]-1,sphereNodes[0,1]+1]).astype(int)

    #Create an empty array for the two end points spaning the user point
	LinePoints=np.zeros((2,2))
	pointPosition=1
    #Store the end points IS THE ORDER IMPORTANT?
	LinePoints[0,:]=contourData[nodeGroup[0],:]
	LinePoints[1,:]=contourData[nodeGroup[pointPosition],:]
	lineCheck=False

	print 'Running lineCheck'
	xPt,yPt = projectPointToLine(LinePoints,point)
	lineCheck=checkPointOnLineInterval(LinePoints,[xPt,yPt])
	if not lineCheck:
		pointPosition=pointPosition+1
		LinePoints[1,:]=contourData[nodeGroup[pointPosition],:]
		xPt,yPt = projectPointToLine(LinePoints,point)
		lineCheck=checkPointOnLineInterval(LinePoints,[xPt,yPt])
		#If the next segment also fales, our picked point is the nearest node
		if not lineCheck:
			xPt=contourData[nodeGroup[0],0]
			yPt=contourData[nodeGroup[0],1]

	#Special case, if the two spanning nodes are the ends nodes
	if ((nodeGroup[0]==0) and (nodeGroup[pointPosition] == len(contourData)-1)) or ((nodeGroup[0]==len(contourData)-1) and (nodeGroup[pointPosition] == 0)):
		print 'Spanning end nodes'
		return [xPt, yPt], len(contourData)-1
	if nodeGroup[0]<nodeGroup[pointPosition]:
		#print 'left is smaller'
		return [xPt, yPt], nodeGroup[0]
	else:
		print 'right is smaller'
		return [xPt, yPt], nodeGroup[pointPosition]

#Function sortContourData takes an input contour data text file, and sorts the file into a list format
def sortContourData(contourData):
	print 'Running sortContourData'
	#loop throught the data set and determine how many contours there are
	
	#Old code block
	#count=0
	#for i in range(len(contourData)-1):
	#	if contourData[i,2] != contourData[i+1,2]:
	#		count=count+1

	#print count
	#print int(contourData[-1,2])
	count = int(contourData[-1,2]-contourData[0,2]+1) #this is silly cause my images start from 0, but my contours start from 1.

	#set up a blank list
	contours=[0]*count
	print '\tNumber of contours:',count
	#set a counter that will count the number of points in a contour
	count=0
	#set a boolean flag toggle to TRUE
	toggle=True
	#zero the contour number data
	contourData[:,2]=contourData[:,2]-contourData[0,2]
	#loop through the contourData list
	for i in range(len(contourData)):
		#set the current working (temporal) contour
		if toggle:
			index=contourData[i,2].astype(int)
			toggle=False

		#print 'i: ',i,'index',contourData[i,2]
		if contourData[i,2]==index:
			count=count+1
		else:
			print 'Contour Number:',index,'Contour Size:',np.shape(contourData[(i-count):i,0:2]),'i:',i,'count',count,'i-count',(i-count)
			contours[index]=contourData[(i-count):i,0:2]
			print contourData[(i-count):i,0:2]
			toggle=True
			count=1

	#you have written some really class A shit code here Stephen, this is gawd aweful, but you 
	#have 4 months left, so we will just role with it, tack on the final index position after
	#loop
	print 'Contour Number:',index,'Contour Size:',np.shape(contourData[(i-count):i,0:2]),'i:',i,'count',count,'i-count',(i-count)
	contours[index]=contourData[(i-count+1):,0:2] #im tired, Im not sure why the plus one is required, basic indexing but meh it works
	print contourData[(i-count+1):,0:2] 
	print contourData[-1,:]
	print contourData[i,:]

	return contours

#Function marches around a contour, placeing the user specified number of slave points at
#equidistance intervals between Master points.
#Dist: the current distance that the crawler is at for a slaveSegment. One it moves to a new slaveSegment, it resest its self to the remainder.
#segDist:The distance of the current dataSegment (not this is a segment between raw data points, not the slave points we are seeding, thats slaveSegment.)
#segSize: The distance of the current slaveSegment, this is some really shitty naming from me, my bad.
#MDist: the current distance for the whole master segment. resets itself when it moves to a new master point
#CPpos: Current control point position (seed points, so could rename this aswell.)
#CNP: Current Node Points, this stores the value of the current node ID from the data array, its used to tell when the crawler has made it to the next masterPoint.
def seedSlavePoints(distVector,contourData,MasterPoints,seedNumber,ax):
	import pdb
	#pdb.set_trace()
	print 'Running seedSlavePoints'
	seedNumber=seedNumber+1
	ControlPoints=np.zeros(((len(MasterPoints)*seedNumber),2))
	print 'Total Number of Control Points:',np.shape(ControlPoints)
	#Start from the first master Point
	CPpos=0
	CNP=MasterPoints[0,2].astype(int)
	#set up a phantom master point again
	MasterPoints=np.concatenate((MasterPoints,MasterPoints[0,:].reshape(1,3)))

	for i in range(len(MasterPoints)-1):
		#set the slaveSegment size for the current slave segment.
		segSize=distVector[i]/(seedNumber)
		#print 'DistVector',distVector[i],'seedNumber',seedNumber-1,'MasterNumber',len(MasterPoints),'segSize',segSize
		#Initislise vaiables and arrays using the current control point (which is the master point for this segment)
		Dist=MDist=nodeCount=0
		#print '\nMpoint:',MasterPoints[i,2],'CPpos:',CPpos,'CNP:',CNP
		ControlPoints[CPpos,:]=MasterPoints[i,:2]
		#increment to the first slave point
		CPpos=CPpos+1
		#looping through data points, while the the current node is not the next master point

		#So I wrote a forwards marching alogrythem CNPi -> CNPi+1 that starts at a point and ends and the same point. this is stupid,
		#because it doesnt work for a single master point, as the first CNP and laster MP are the same, got 4 months left,
		#dont wont to break all the code by trying to fix this, so im adding in a nasty fix to allow a special case for
		#single MP. This whole thing can be fixed by changing it to a system that sets CNP to MP(1)+1 and then does CNP<-CNP+1
		#(looks backwards).

		oneMPcase=False
		if len(MasterPoints)==2: #so if there is only 1 MP and then the ghost point
			oneMPcase=True

		while (CNP!=MasterPoints[i+1,2] or oneMPcase):
			oneMPcase=False #switch this stupid thing off now
			#If else to check if we are at the end of the array, to loop back to the beginning.
			if CNP==(len(contourData)-1):
				target=0
			else:
				target=CNP+1

			#calculate the length of the data segment that we are on.
			segDist=calculateLineLength(contourData[CNP,:],contourData[target,:])
			MDist=MDist+segDist
			#print '\tCPpos',CPpos,'CNP',CNP,'SegDist',segDist,'Dist',Dist,'MDist',MDist

			#if the current distance, plus the length of the data segment is larger then our current slaveSegment spacing. we will calculate where the slave seed should be placed, and what remains
			if Dist+segDist >= segSize:
				#print 'Mpoint:',MasterPoints[i,2],'CPpos:',CPpos,'CNP:',CNP
				ControlPoints[CPpos,:],remainder = calculatePointPositionOnLine(segSize-Dist,segDist,contourData[CNP,:],contourData[target,:])
				nodeCount=nodeCount+1
				#print '\t\tplacing a seed point at position',CPpos,'Point Positions [',ControlPoints[CPpos,0],ControlPoints[CPpos,1],'], remainder:',remainder
				Dist=remainder
				CPpos=CPpos+1
				#on the offchance that the segment distance of the data points is huge, begger then our seed segment size, we nee to place multiple points on the data segment.
				while remainder > segSize:
					if nodeCount==seedNumber-1:
						remainder=0
					else:
						print 'running while loop: Remainder:', remainder, 'segSize:',segSize
						print 'Mpoint:',MasterPoints[i,2],'CPpos:',CPpos,'CNP:',CNP
						ControlPoints[CPpos,:],remainder = calculatePointPositionOnLine((segDist-(remainder-segSize)),segDist,contourData[CNP,:],contourData[target,:])
						nodeCount=nodeCount+1
						CPpos=CPpos+1

				Dist=remainder
			else:
				Dist=Dist+segDist

			if nodeCount==seedNumber-1:
				CNP=MasterPoints[i+1,2].astype(int)
			else:
				CNP=target
			#print 'CNP value:',CNP,'Distance',Dist

	return ControlPoints

def calculatePointPositionOnLine(segSize,segDist,point1,point2):
	#calculate the fraction along the line length that the seed point will lie
	#frac=((segSize-segDist)/segSize)
	frac=(segSize/segDist)
	xPos=(point2[0]-point1[0])*frac + point1[0]
	yPos=(point2[1]-point1[1])*frac + point1[1]
	remainder=(1-frac)*segDist
	#print '\t\tRunning calculatePointPositionOnLine'
	#print '\t\tFrac value:',frac,'remainder ratio',(1-frac)
	return [xPos,yPos],remainder


#Function takes a data intervale of 2D coordinates, and marches along calculating the total distance.
def calculateDistanceVector(contourData,MasterPoints):
	print 'Running calculateDistanceVector'
	distVector=np.zeros(len(MasterPoints))
	#Start from the first master Point, check that its not right at the end
	#print 'Contour Data:\n', contourData
	if MasterPoints[0,2]==(len(contourData)-1):
		CNP=0
	else:
		CNP=MasterPoints[0,2].astype(int)+1

	MPcount=0
	distVector[MPcount]=calculateLineLength(contourData[MasterPoints[0,2].astype(int),:],contourData[CNP,:])
	#create a phantom node to stop the loop condition
	MasterPoints=np.concatenate((MasterPoints,MasterPoints[0,:].reshape(1,3)))
	while CNP!=MasterPoints[0,2].astype(int):
		#check for a loop around
		#print 'Max/CNP:',(len(contourData)-1),'/',CNP
		if CNP==(len(contourData)-1):
			distVector[MPcount]=distVector[MPcount]+calculateLineLength(contourData[CNP,:],contourData[0,:])
			#print 'nodes (',CNP,') [',contourData[CNP,0],contourData[CNP,1],'] (',0,') [',contourData[0,0],contourData[0,1],'] Distance',distVector[MPcount]
			CNP=0
		else:
			distVector[MPcount]=distVector[MPcount]+calculateLineLength(contourData[CNP,:],contourData[CNP+1,:])
			#print 'nodes (',CNP,') [',contourData[CNP,0],contourData[CNP,1],'] (',CNP+1,') [',contourData[CNP+1,0],contourData[CNP+1,1],'] Distance',distVector[MPcount]
			CNP=CNP+1

		#print 'Distance',distVector[MPcount]
		if (CNP==MasterPoints[MPcount+1,2]):
			MPcount=MPcount+1
			#print 'MPcount',MPcount

	return distVector



def calculateLineLength(Point1,Point2):
	return np.sqrt((Point2[0]-Point1[0])**2 + (Point2[1]-Point1[1])**2)

def projectPointToLine(Line,Point):
	print 'Running projectPointToLine'
	#print 'Line:',Line,'Point1:',Point[0],'Point2',Point[1]
	#Check if the line is vertical
	if abs(Line[1,0]-Line[0,0])<0.0001:
		print 'Vertical line detected'
		xPt=Line[1,0]
		yPt=Point[1]
	elif abs(Line[1,1]-Line[0,1])<0.0001:
		print 'Horizontal line detected'
		xPt=Point[0]
		yPt=Line[1,1]
	else:
		#Calculate the gradient Mt (target) and constant Ct
		Ct=(Line[0,1]*Line[1,0] - Line[1,1]*Line[0,0])/(Line[1,0]-Line[0,0])
		Mt=(Line[0,1]-Ct)/Line[0,0]
		#The gradient of the orthogonal line is the inverse of Mt
		Ms=-1/Mt
		#Given the point and the gradient, calculate the constant
		Cs=Point[1]-Point[0]*Ms
		#Calculate the intercept between the two lines
		xPt=(Cs-Ct)/(Mt-Ms)
		yPt=Mt*xPt+Ct

	return xPt,yPt

def checkPointOnLineInterval(Line,Point):
	yChecks=[Line[0,1]-Point[1],Line[1,1]-Point[1]]
	xChecks=[Line[0,0]-Point[0],Line[1,0]-Point[0]]
	passY=False
	passX=False
	if (np.sign(yChecks[0])!=np.sign(yChecks[1])):
		passY=True
	if (np.sign(xChecks[0])!=np.sign(xChecks[1])):
		passX=True

	#Return results
	if (passY and passX):
		return True
	else:
		return False

def readContourFile(filePath):
	#Function readContourFile is used to read in createControlPoints.py write files.
	#These files contain two bits of header information. Their Z position, and 
	#the number of contours in the file. These are first read and then np.loadtxt 
	#and reshpare are used to format the contour data.
	print 'Running readContourFile...'
	print 'Opening from location',filePath
	f=open(filePath,'r')
	#get the Z position
	zPos=f.readline().split()[2]
	#print 'zPos:',zPos
	#get the number of contours position
	numContours=f.readline().split()[2]
	#print 'Number of Contours:',numContours

	#read in the contour data, reshape it and store in a list
	contours=np.loadtxt(filePath)
	#print 'reshaping',contours.shape,'to [',int(numContours),contours.shape[0]/int(numContours),2,']'
	contours=contours.reshape((int(numContours),contours.shape[0]/int(numContours),2))
	contours.tolist()

	return zPos,numContours,contours

def createContractionMatrix(activeContours,contourInformation,contourData,CMap):

	print 'Running createContractionMatrix...'
	
	#Print out some sanity statistics
	print 'Contour information:'
	for i in range(np.size(activeContours)):
		print '\t',activeContours[i].rstrip(),'zPos:',contourInformation[i][0],'number Contours',contourInformation[i][1]

	#So the current plan is, to bit a time range (roughly 2 minutes).
	#we create a nxd array where n is the number of slice planes we are using, and d is seconds.
	#using our contraciton map from literature, we will place active contours in the correct
	#times, and duplicate their first and last slices to rill out time on either side.
	#will need to calculate a suitable end time, as some contours contract for a while.
	#This will be dependent on when they start in the contraciton map. 
	#For now we will take the largest contraciton time and tack it onto the end time in the 
	#CMAP
	endtime=np.int(CMap[-1,1]+max(np.asarray(contourInformation)[:,1].astype(float)))+1
	print 'End time is:',endtime

	#set up the contractionMatrix
	#contractionMatrix=[[0]*np.shape(activeContours)[0]]*endtime <-never do this
	contractionMatrix=[[0]*np.shape(activeContours)[0] for i in range(endtime)]
	#
	print 'contractionMatrix shape:',np.shape(contractionMatrix)
	MapFunc=interp1d(CMap[:,0],CMap[:,1],kind='cubic')
	print CMap

	for i in range(np.shape(activeContours)[0]):
		print np.float(contourInformation[i][0])
		spoint=np.int(MapFunc(np.float(contourInformation[i][0])))
		print 'Processing slice:',i,'at',contourInformation[i][0]
		print '\tStart Points:',spoint, 'End Point:',spoint+int(contourInformation[i][1])
		print 'looping from',0,'to',spoint
		for j in range(spoint):
			#print 'i:',i,'j:',j,'slice:',0
			contractionMatrix[j][i]=contourData[i][0]

		print 'looping from',spoint,'to',np.int(spoint+int(contourInformation[i][1]))
		for j in range(spoint,np.int(spoint+int(contourInformation[i][1]))):
			#print 'i:',i,'j:',j,'slice',(j-spoint)
			contractionMatrix[j][i]=contourData[i][j-spoint]

		print 'looping from',np.int(spoint+int(contourInformation[i][1])),'to',endtime
		for j in range(np.int(spoint+int(contourInformation[i][1])),endtime):
			#print 'i:',i,'j:',j,'slice: end'
			contractionMatrix[j][i]=contourData[i][-1]

	return contractionMatrix,endtime

def calculateCentroidPosition(x,y):
	return sum(x)/len(x),sum(y)/len(y)

def calculateCentroidSize(x,y):
	cX,cY=calculateCentroidPosition(x,y)
	#print cX
	#print cY
	#print sum((x-cX)**2 + (y-cY)**2)
	#print np.sqrt(sum((x-cX)**2 + (y-cY)**2))
	return np.sqrt(sum((x-cX)**2 + (y-cY)**2))

def calculateClosedContourArea(x,y):
	return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

def PolygonArea(x,y):
    n = len(x) # of corners
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += x[i] * y[j]
        area -= x[j] * y[i]
    area = abs(area) / 2.0
    return area

    #https://stackoverflow.com/questions/24467972/calculate-area-of-polygon-given-x-y-coordinates <-test these later dawg

def basicCubicSpline(n,xn,yn):
	#Function basicCubicSpline takes the number of data points (n), and their x and y values, and computes 
	#equations for basic cubic splines, where the data points are the knots
	#Coefficents are returned in the form d*x^3 + c*x^2 + b*x + yn

	#set up a series of constant vectors
	h = np.zeros(n-1) #x spacing between knots
	alpha = np.zeros(n-1)
	l = np.zeros(n+1)
	u = np.zeros(n)
	z = np.zeros(n+1)
	b = np.zeros(n)
	c = np.zeros(n+1)
	d = np.zeros(n)

	#Calculate X spacing
	for i in range(n-1):
		h[i] = xn[i+1]-xn[i]

	#Calculate alpha values
	for i in range(1,n-1):
		alpha[i] = (3./h[i])*(yn[i+1]-yn[i])-(3./h[i-1])*(yn[i]-yn[i-1])

	#Set end point conditions for a natural spline
	l[0] = l[n] = 1
	u[0] = c[n] = 0
	z[0] = z[n] = 0 

	#Calculate constants l,u,z
	for i in range(1,n-1):
		l[i] = 2*(xn[i+1]-xn[i-1])-h[i-1]*u[i-1]
		u[i] = h[i]/l[i]
		z[i] = (alpha[i]-h[i-1]*z[i-1])/l[i]

	for i in range(n-2,-1,-1):
		c[i] = z[i] - u[i]*c[i+1]
		b[i] = (yn[i+1]-yn[i])/h[i]-h[i]*(c[i+1]+2*c[i])/3
		d[i] = (c[i+1]-c[i])/(3*h[i])

	return yn,b,c,d



def basicCubicSplinePlot(a,b,c,d,xn):
	#Function basicCubicSplinePlot is used to plot spines from function basicCubicSpline for test purposes. It takes 
	#the equation coefficents in the form d*x^3 + c*x^3 + b*x + a and the xn positions
	n = np.size(xn)
	if (np.size(a) != n) or (np.size(b) != n) or (np.size(c)-1 != n) or (np.size(d) != n): 
		print 'Error: vector lengths do not match'
		return 0
	else:
		for i in range(n-1):
			root = poly1d(xn[i],True)
			poly = 0
			poly = d[i]*(root)**3
			poly = poly + c[i]*(root)**2
			poly = poly + b[i]*root
			poly = poly + a[i]
			#print 'poly:',poly
			#set up a vector of x values from xi to xi+1
			x_step=0.05
			x_pts = np.arange(xn[i],xn[i+1]+x_step,x_step)
			plt.plot(x_pts,poly(x_pts),'-')
		return 1

def transformCubicCoeffiencets(a,b,c,d,xn):
	#function uses poly1d to return the coefficents of a spline in the form of a cubic equation.
	n = np.size(xn)
	aCoeffs=np.zeros(n)
	bCoeffs=np.zeros(n)
	cCoeffs=np.zeros(n)
	dCoeffs=np.zeros(n)
	if (np.size(a) != n) or (np.size(b) != n) or (np.size(c)-1 != n) or (np.size(d) != n): 
		print 'Error: vector lengths do not match'
		return 0,0,0,0
	else:
		for i in range(n-1):
			root = poly1d(xn[i],True)
			poly = 0
			poly = d[i].tolist()*(root)**3
			poly = poly + c[i].tolist()*(root)**2
			poly = poly + b[i].tolist()*root
			poly = poly + a[i].tolist()

			tempCoeffs=np.zeros(4)
			tempCoeffs[(4-(poly.order+1)):]=poly.c

			aCoeffs[i],bCoeffs[i],cCoeffs[i],dCoeffs[i]=tempCoeffs
			#aCoeffs[i],bCoeffs[i],cCoeffs[i],dCoeffs[i]=poly.c
		return aCoeffs,bCoeffs,cCoeffs,dCoeffs






