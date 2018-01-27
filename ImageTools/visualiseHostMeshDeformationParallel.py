#script is used to visualise the host mesh as it deforms over time.
from traits.api import HasTraits, Int, Range, Bool, Button, Str, Enum, List, Instance, on_trait_change
from traitsui.api import View, Item, HGroup, HSplit, VGroup, Tabbed, EnumEditor, TextEditor, Group
from mayavi.core.pipeline_base import PipelineBase
from mayavi.core.ui.api import SceneEditor, MlabSceneModel, MayaviScene
import numpy as np

print 'RUnning visualiseHostMeshDeformationParallel...'

workingDir='/home/stephen/OpenFOAM/Simulations2/Meshing/SuperRoofMeshes/SuperRoofSmoothedReticulum'
patchName='wall'
numProc=10 

#pointTemp=np.loadtxt(workingDir+'/constant/patchPositions/' + patchName + '/' + str(0)  + '/patchDisplacements')[::10]*1000
pointTemp=np.zeros((1,3))
for i in range(numProc):
	#read in the patch file from each processor and smoosh them togeather
	temp=np.loadtxt(workingDir+'/processor'+str(i)+'/constant/patchPositions/' + patchName + '/' + str(0)  + '/patchDisplacements')
	pointTemp=np.vstack([pointTemp,np.asarray(temp)])

pointTemp=np.delete(pointTemp,(0),0)[::20]*1000
print np.shape(pointTemp)


sourceTemp=np.loadtxt(workingDir+'/HMresults/'+ str(1) + '/source_points_fitting_hmf')
targetTemp=np.loadtxt(workingDir+'/HMresults/'+ str(1) + '/target_points')
HMOTemp=np.loadtxt(workingDir+'/HMresults/'+ str(1) + '/HostMeshOrig')
HMDTemp=np.loadtxt(workingDir+'/HMresults/'+ str(1) + '/HostMeshDeform')

numT=8

pointData=np.zeros((pointTemp.shape[0],pointTemp.shape[1],numT))
sourceData=np.zeros((sourceTemp.shape[0],sourceTemp.shape[1],numT))
targetData=np.zeros((targetTemp.shape[0],targetTemp.shape[1],numT))
HMOData=np.zeros((HMOTemp.shape[0],HMOTemp.shape[1],numT))
HMDData=np.zeros((HMDTemp.shape[0],HMDTemp.shape[1],numT))
for i in range(numT):
	print 'loading time:', i
	#pointData[:,:,i]=np.loadtxt(workingDir+'/constant/patchPositions/' + patchName + '/' + str(i+1)  + '/patchDisplacements')[::10]*1000
	sourceData[:,:,i]=np.loadtxt(workingDir+'/HMresults/'+ str(i+1) + '/source_points_fitting_hmf')
	targetData[:,:,i]=np.loadtxt(workingDir+'/HMresults/'+ str(i+1) + '/target_points')
	HMOData[:,:,i]=np.loadtxt(workingDir+'/HMresults/'+ str(i+1) + '/HostMeshOrig')
	HMDData[:,:,i]=np.loadtxt(workingDir+'/HMresults/'+ str(i+1) + '/HostMeshDeform')

	#now we do the pointNodes
	pointTemp=np.zeros((1,3))
	for t in range(numProc):
		#read in the patch file from each processor and smoosh them togeather
		temp=np.loadtxt(workingDir+'/processor'+str(t)+'/constant/patchPositions/' + patchName + '/' + str(i)  + '/patchDisplacements')
		pointTemp=np.vstack([pointTemp,np.asarray(temp)])

	print 'check'
	print np.shape(np.delete(pointTemp,(0),0))
	print np.shape(pointData)
	FUCK=np.delete(pointTemp,(0),0)
	pointData[:,:,i]=FUCK[::20]*1000


class TestModel(HasTraits):
	scene = Instance(MlabSceneModel, ())
	plot = Instance(PipelineBase)
	print 'life is pain'
	t=int(0)
	s=1000

	#self.scene.mlab.points3D(testData[:,0,t]*s,testData[:,1,t]*s,testData[:,2,t]*s,'ob')
	@on_trait_change('t')
	def update_plot(self):
		if self.plot is None:
			print 'FIRST PLOT'
			self.plot=self.scene.mlab.points3D(testData[:,0,t]*s,testData[:,1,t]*s,testData[:,2,t]*s,'ob')
		else:
			print 'T has changed to:',t
			self.plot=self.scene.mlab.points3D(testData[:,0,t]*s,testData[:,1,t]*s,testData[:,2,t]*s,'rb')
	
	#@on_trait_change('t')
	#def update_plot(self):
	#	if self.plot is None:
	#		print 'FIRST PLOT'
	#		self.plot=self.scene.mlab.points3D(testData[:,0,t]*s,testData[:,1,t]*s,testData[:,2,t]*s,'ob')
	#	else:
	#		print 'T has changed to:',t
	#		self.plot=self.scene.mlab.points3D(testData[:,0,t]*s,testData[:,1,t]*s,testData[:,2,t]*s,'rb')

	view=View(Item('scene',editor=SceneEditor(scene_class=MayaviScene),height=1000,width=2000,show_label=False),Group('_','t'),resizable=True,)

#Test=TestModel()
#Test.configure_traits()

class testVI(HasTraits):

	renderAll = Button()
	renderImagePlane = Button()
	renderData = Button()
	
	imageList0 = 'None'
	imageList = List([])
	imagePlane = Enum( ('x_axes', 'y_axes', 'z_axes') )
	imageVisible = Bool(False)
	
	dataList0 = 'None'
	dataList = List([])
	dataVisible = Bool(True)
	dataUpdate = Button()
		
	saveImageFilename = Str('screenshot.jpeg')
	saveImageWidth = Str('2000')
	saveImageLength = Str('1500')
	saveImage = Button()
	
	
	scene = Instance(MlabSceneModel, ())

	t=Int(3)
	plot = Instance(PipelineBase)

	colours = {'bone':(0.84705882, 0.8, 0.49803922)}
	defaultColor = colours['bone']	# bone
	#~ defaultColor = (0.5,0.5,0.5)	# bone

	view = View(HSplit(
                    Tabbed(
	                    VGroup(
	                        Item('renderAll', show_label=False, label='render all objects', springy=True),
	                        #~ Item('renderImagePlane', show_label=False, label='render image plane', springy=True),
	                       
	                        Item('renderData', show_label=False, label='render data clouds', springy=True),
	                        
	                        VGroup(	
								Item( 'imageList0', show_label=True, label='images', editor = EnumEditor(name='imageList') ),
								HGroup(
									Item( 'imageVisible', show_label=True, label='visible', springy=True ),
									Item( 'imagePlane', show_label=True, label='image plane', springy=True),
									)
								),
								
	                        VGroup(
								Item( 'dataList0', show_label=True, label='Dataclouds',editor = EnumEditor(name='dataList') ),
								HGroup(
									Item( 'dataVisible', show_label=True, label='visible', springy=True ),
									Item( 'dataUpdate', show_label=False, label='update', springy=True )
									)
								),
								
							VGroup(
								Item( 'saveImageFilename', show_label=True, label='filename', editor = TextEditor() ),
								HGroup(
									Item( 'saveImageWidth', show_label=True, label='W', editor = TextEditor() ),
									Item( 'saveImageLength', show_label=True, label='L', editor = TextEditor() ),
									),
								Item( 'saveImage', show_label=False, label='saveImage', springy=True ),
								),
								
	                        label='Render'
	                    ),
	                                        
                    springy=False),
                    
                    Item('scene', editor=SceneEditor(scene_class=MayaviScene), height=600, width=800, show_label=False),
                ),
                resizable=True,
                )

	#@on_trait_change('t')
	#def _drawData( self, name ):
	#	self.scene.disable_render = True
	#	d = self.data[name]
	#	s = self.dataScalar.get(name)
	#	renderArgs = self.dataRenderArgs[name]
	#	if s!=None:
			#~ self.sceneObjectData[name] = self.scene.mlab.points3d( d[:,0], d[:,1], d[:,2], s,
															   #~ mode='point', scale_factor=0.5,
															   #~ scale_mode='none',
															    #~ name=name )
	#		self.sceneObjectData[name] = self.scene.mlab.points3d( d[:,0,self.t], d[:,1,self.t], d[:,2,self.t], s, **renderArgs)
	#	else:
			#~ self.sceneObjectData[name] = self.scene.mlab.points3d( d[:,0], d[:,1], d[:,2],
															   #~ mode='point', scale_factor=0.5,
															   #~ color=(0.0,1.0,0.0),
															    #~ name=name )
	#		self.sceneObjectData[name] = self.scene.mlab.points3d( d[:,0,self.t], d[:,1,self.t], d[:,2,self.t], **renderArgs )
		
	#	self.scene.disable_render = False


	def __init__( self ):
		HasTraits.__init__(self)
		
		self.I = None
		
		self.images = {}
		self.imageCounter = 0
		self.sceneObjectImages = {}
		self.imageRenderArgs = {}
		
		self.data = {}
		self.dataScalar = {}
		self.dataCounter = 0
		self.sceneObjectData = {}
		self.dataRenderArgs = {}
		
		
		self.sceneObjectImageVolume = None
		self.modeVectorFlip = 1.0
			
		self.view = None
		self.roll = None

		self.onCloseCallback = None

		# image plane widget attributes
		self._ipw_picked_obj = None
		self._ipw_picked_points = []

	def addOnCloseCallback(self, callback):
		self.onCloseCallback = callback

	def closed(self, info, is_ok):
		if self.onCloseCallback!=None:
			self.onCloseCallback()

	def addData( self, name, d, scalar=None, renderArgs=None ):
		
		if name not in self.data.keys():
			self.dataCounter += 1
			self.dataList.append( name )
						
		self.data[name] = d
		if renderArgs==None:
			self.dataRenderArgs[name] = {}
		else:
			self.dataRenderArgs[name] = renderArgs
		if scalar is not None:
			self.dataScalar[name] = scalar
		else:
			self.dataScalar[name] = None
			
		self.dataList0 = name

	def _drawData( self, name ):
		self.scene.disable_render = True
		d = self.data[name]
		s = self.dataScalar.get(name)
		t = self.t
		renderArgs = self.dataRenderArgs[name]
		if s!=None:
			#~ self.sceneObjectData[name] = self.scene.mlab.points3d( d[:,0], d[:,1], d[:,2], s,
															   #~ mode='point', scale_factor=0.5,
															   #~ scale_mode='none',
															    #~ name=name )
			self.sceneObjectData[name] = self.scene.mlab.points3d( d[:,0], d[:,1], d[:,2], s, **renderArgs)
		else:
			#~ self.sceneObjectData[name] = self.scene.mlab.points3d( d[:,0], d[:,1], d[:,2],
															   #~ mode='point', scale_factor=0.5,
															   #~ color=(0.0,1.0,0.0),
															    #~ name=name )
			self.sceneObjectData[name] = self.scene.mlab.points3d( d[:,0], d[:,1], d[:,2], **renderArgs )
		
		self.scene.disable_render = False

	def _dataVisible_changed( self ):
		try:
			self.sceneObjectData[self.dataList0].visible = self.dataVisible
		except KeyError:
			if self.dataVisible:
				self._drawData( self.dataList0 )

	def _dataUpdate_fired( self ):
		self.updateData( self.dataList0 )
		
	def updateData( self, name, coords=None ):
		if coords is None:
			d = self.data[name]
		else:
			d = coords
		s = self.dataScalar.get(name)
		renderArgs = self.dataRenderArgs[name]
		t=self.t
		
		try:
			if s is not None:
				self.sceneObjectData[name].actor.mapper.scalar_visibility=True
				self.sceneObjectData[name].mlab_source.reset( x=d[:,0], y=d[:,1], z=d[:,2], s=s )
			else:
				if 'color' not in renderArgs:
					color = self.defaultColor
				else:
					color = renderArgs['color']
					
				self.sceneObjectData[name].actor.mapper.scalar_visibility=False
				self.sceneObjectData[name].actor.property.specular_color = color
				self.sceneObjectData[name].actor.property.diffuse_color = color
				self.sceneObjectData[name].actor.property.ambient_color = color
				self.sceneObjectData[name].actor.property.color = color
				self.sceneObjectData[name].mlab_source.reset( x=d[:,0], y=d[:,1], z=d[:,2] )
		except KeyError:
			self._drawData( name )
			
	#def _renderData_fired( self ):
		""" redraw all data clouds
		"""
	#	for i in self.sceneObjectData.values():
	#		i.remove()
			
	#	self.sceneObjectData = {}
	#	for i in self.data.keys():
	#		self._drawData( i )



v=testVI()
v.addData('Origin',np.zeros((4,3)),renderArgs={'mode':'sphere', 'scale_factor':5, 'color':(1,0,0)})
#v.addData('my_array',my_array,renderArgs={'mode':'sphere', 'scale_factor':3, 'color':(0.5,0.5,0.5)})
for i in range(numT):
	v.addData('pointData'+str(i+1),pointData[:,:,i],renderArgs={'mode':'sphere', 'scale_factor':2, 'color':(0.1,0.1,0.1)})

for i in range(numT):	
	v.addData('sourceData'+str(i+1),sourceData[:,:,i],renderArgs={'mode':'sphere', 'scale_factor':3, 'color':(0.8,0.1,0.1)})

for i in range(numT):	
	v.addData('targetData'+str(i+1),targetData[:,:,i],renderArgs={'mode':'sphere', 'scale_factor':3, 'color':(0.1,0.1,0.8)})

for i in range(numT):	
	v.addData('HM0'+str(i+1),HMOData[:,:,i],renderArgs={'mode':'sphere', 'scale_factor':3, 'color':(0.5,0.1,0.1)})

for i in range(numT):	
	v.addData('HMD'+str(i+1),HMDData[:,:,i],renderArgs={'mode':'sphere', 'scale_factor':3, 'color':(0.1,0.1,0.5)})

#v.addData('testData'+str(1),pointData[:,:,0],renderArgs={'mode':'sphere', 'scale_factor':1, 'color':(0.1,0.1,0.1)})
#v.addData('testData'+str(2),pointData[:,:,1],renderArgs={'mode':'sphere', 'scale_factor':1, 'color':(0.1,0.1,0.1)})
#v.addData('testData'+str(3),pointData[:,:,2],renderArgs={'mode':'sphere', 'scale_factor':1, 'color':(0.1,0.1,0.1)})
#v.addData('testData'+str(4),pointData[:,:,3],renderArgs={'mode':'sphere', 'scale_factor':1, 'color':(0.1,0.1,0.1)})


v.configure_traits()
v.scene.background=(0,0,0)


