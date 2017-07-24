#List of function tools for the sementor

import sys
import numpy as np

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





					




