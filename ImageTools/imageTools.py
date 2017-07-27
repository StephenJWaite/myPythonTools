#List of function tools for the sementor

import sys
import numpy as np
import matplotlib.pyplot as plt

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
			print 'Point ', i, ': ', contourData[i,:]
			dist=np.sqrt((contourData[i,0]-point[0])**2 + (contourData[i,1]-point[1])**2)
			print 'Dist: ',dist
			if sphereNodes[0,0]==-1:
				sphereNodes[0,:]=[dist,i]
			else:
				sphereNodes=np.append(sphereNodes,[[dist,i]],axis=0)

	#sort the sphereNodes
	print 'Dist vector shape: ', np.shape(sphereNodes)
	print 'Dist vector:', sphereNodes
	sphereNodes=sphereNodes[np.argsort(sphereNodes[:,0]),:]
	print 'Sorted Dist vector:', sphereNodes
	#Calculate the distance between the two closets nodes
	#nodeDist=(np.sqrt((contourData[sphereNodes[0,1],0]
	#          -contourData[sphereNodes[1,1],0])**2
	#          +(contourData[sphereNodes[0,1],1]
	#          -contourData[sphereNodes[1,1],1])**2
	#         ))
	#Calculate the angle between contour nodes and the picked node
	LinePoints=np.zeros((2,2))
	print 'LinePoints:', LinePoints
	print 'sphereNodes:', sphereNodes[0,1]
	print 'contourData[sphereNodes[0,1],:]:',contourData[sphereNodes[0,1].astype(int),:]
	LinePoints[0,:]=contourData[sphereNodes[0,1].astype(int),:]
	LinePoints[1,:]=contourData[sphereNodes[1,1].astype(int),:]
	xPt,yPt = projectPointToLine(LinePoints,point)
	print 'xPt: ',xPt, ' yPlst: ',yPt


def projectPointToLine(Line,Point):
	print 'Running projectPointToLine'
	print 'Line shape: ', np.shape(Line)
	print 'Line values: ', Line
	print 'Point values: ', Point
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
	print 'xPt: ',xPt, ' yPt: ',yPt

	#linePlot = plt.figure()
	#axLine = linePlot.add_subplot(111)
	#axLine.plot([Line[0,0],Line[1,0]],[Line[0,1],Line[1,1]])
	#axLine.scatter(Point[0],Point[1], facecolors='g', edgecolors='g')
	#axLine.scatter(xPt,yPt, facecolors='r', edgecolors='r')
	#axLine.plot([Point[0],xPt],[Point[1],yPt],color='r')
	#axLine.axis('equal')

	return xPt,yPt









					




