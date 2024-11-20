#!/usr/bin/env python
import sys
sys.path.append( '../stochastic')
import vtk
import numpy as np
from scipy import interpolate
from vtk.util.colors import *
import random
math = vtk.vtkMath()


class GraphicObject(object):
	pass

class PanelViewer3D(object):

	def __init__(self):
		self.output=GraphicObject()
		self.output.quadCount=-1
		self.output.quadl=[]
		self.output.quads=vtk.vtkCellArray()
		self.output.points = vtk.vtkPoints();
		self.output.polyVerticesCount=0
		self.output.polyVertices= vtk.vtkPoints();
		self.polyCount=0
		self.output.polyList=[]
		self.output.polygons=vtk.vtkCellArray()
		self.output.objects=vtk.vtkAppendPolyData()
		
		self.axisSystems= []
		self.pointSets = []
		self.polylines=[]
		self.lines=[]
		self.textLabels=[]
		self.spheres=[]
		self.arrows=[]
		self.textScale=5e-3
		self.colors=vtk.vtkFloatArray()
		self.scG=vtk.vtkFloatArray()
		self.scGid=vtk.vtkFloatArray()
		self.scPGid=vtk.vtkFloatArray()
		self.force=vtk.vtkFloatArray()
		self.force.SetName('Force')
		self.force.SetNumberOfComponents(3)
		
		self.scS=vtk.vtkFloatArray()
		self.colors.SetName('CP')
		self.scG.SetName('G')
		self.scGid.SetName('gid')
		self.scPGid.SetName('PgId')

		self.scS.SetName('S')
		self.scalingFactor=1.0
		self.showAxisSystems=True#False
		self.showCompass=True
		self.showAxes=True
		self.showScalarBar=True
		self.showCaption=True
		self.showXYPlane=True
		self.showYZPlane=True
		self.showZXPlane=True
		self.xlabel='x'
		self.ylabel='y'
		self.zlabel='z'
		self.xrange=None
		self.yrange=None
		self.zrange=None
		self.nbSize=80
		self.sampleSpacing=8.0
		self.caption=''
		self.camera = vtk.vtkCamera()
		
		self.writePolydata=False
		self.writeObjects=False

		self.ren = vtk.vtkRenderer();
		self.renWin = vtk.vtkRenderWindow();
		self.renWin.AddRenderer(self.ren);
		self.renWin.SetWindowName("VLM-viewer v0.1 ")
		self.iren = vtk.vtkRenderWindowInteractor();
		self.iren.SetRenderWindow(self.renWin);
		style1 = vtk.vtkInteractorStyleTrackballCamera()
		self.iren.SetInteractorStyle(style1)
		self.ren.SetBackground(1, 1, 1);
		self.renWin.SetSize(800, 600);

	def show(self):
		self.ren.ResetCamera()
		self.iren.Initialize();
		self.renWin.Render();
		self.iren.Start();
		self.close()


	def close(self):
		self.renWin.Finalize()
		self.iren.TerminateApp(); 
		del self.ren
		del self.renWin 
		del self.iren

	def setCaption(self,text):
		self.caption=text
		self.renWin.SetWindowName(" " +str(text))

	def setScalingFactor(self,sf):
		self.scalingFactor=sf

	def setNeighborhoodSize(self,nbSize):
		self.nbSize=nbSize

	def setSampleSpacing(self,sampleSpacing):
		self.sampleSpacing=sampleSpacing

	def writeVTP(self,fname,i=0):
		writer=vtk.vtkXMLPolyDataWriter()
		writer.SetFileName(fname)
		
		appendFilter=vtk.vtkAppendPolyData()
		if self.writePolydata:
			appendFilter.AddInputData(self.output.polydata)
		if self.writeObjects:
			appendFilter.AddInputConnection(self.output.objects.GetOutputPort())
		writer.SetInputConnection(appendFilter.GetOutputPort())
		writer.Write()
	
	def writeVTK(self,fname,i=0):
		writer=vtk.vtkPolyDataWriter()
		writer.SetFileName(fname)
		
		appendFilter=vtk.vtkAppendPolyData()
		if self.writePolydata:
			appendFilter.AddInputData(self.output.polydata)
		if self.writeObjects:
			appendFilter.AddInputConnection(self.output.objects.GetOutputPort())
		writer.SetInputConnection(appendFilter.GetOutputPort())
		writer.Write()

	def newPolyLine(self):
		self.pointSets.append( vtk.vtkPoints())
		self.polylines.append( vtk.vtkPolyLine())
		return len(self.pointSets) -1

	def addTextLabel(self,xyz,txt):
		l={}
		l['xyz']=xyz
		l['txt']=txt
		self.textLabels.append(l)


	def addLine(self,P1,P2):
		line = vtk.vtkLineSource()
		line.SetPoint1(P1[0],P1[1],P1[2])
		line.SetPoint2(P2[0],P2[1],P2[2])

		mapper = vtk.vtkPolyDataMapper()
		source=line
		if vtk.VTK_MAJOR_VERSION <= 5:
			mapper.SetInput(source.GetOutput())
		else:
			#mapper.SetInputData(source)
			mapper.SetInputConnection(source.GetOutputPort())
		
		actor = vtk.vtkActor()
		actor.SetMapper(mapper)

		actor.GetProperty().SetColor(0,0,0)
		self.lines.append(actor)
		self.output.objects.AddInputConnection(source.GetOutputPort())
		
	
	def addArrow(self,P1,P2,scale=1.0,c=(1.,.3,.3)):
		arrowSource = vtk.vtkArrowSource()
		arrowSource.SetShaftRadius(0.01)
		arrowSource.SetTipRadius(0.03)
		arrowSource.SetTipLength(0.2)
		startPoint=P1; endPoint=P2;
		# Compute a basis
		normalizedX = [0 for i in range(3)]
		normalizedY = [0 for i in range(3)]
		normalizedZ = [0 for i in range(3)]
		
		# The X axis is a vector from start to end
		math.Subtract(endPoint, startPoint, normalizedX)
		length = math.Norm(normalizedX)
		math.Normalize(normalizedX)
		
		# The Z axis is an arbitrary vector cross X
		arbitrary = [0 for i in range(3)]
		arbitrary[0] = random.uniform(-10,10)
		arbitrary[1] = random.uniform(-10,10)
		arbitrary[2] = random.uniform(-10,10)
		math.Cross(normalizedX, arbitrary, normalizedZ)
		math.Normalize(normalizedZ)
		
		# The Y axis is Z cross X
		math.Cross(normalizedZ, normalizedX, normalizedY)
		matrix = vtk.vtkMatrix4x4()
		
		# Create the direction cosine matrix
		matrix.Identity()
		for i in range(3):
			matrix.SetElement(i, 0, normalizedX[i])
			matrix.SetElement(i, 1, normalizedY[i])
			matrix.SetElement(i, 2, normalizedZ[i])
		
		# Apply the transforms
		transform = vtk.vtkTransform()
		transform.Translate(startPoint)
		transform.Concatenate(matrix)
		transform.Scale(length*scale, length*scale, length*scale)
		
		# Transform the polydata
		transformPD = vtk.vtkTransformPolyDataFilter()
		transformPD.SetTransform(transform)
		transformPD.SetInputConnection(arrowSource.GetOutputPort())
		
		self.output.objects.AddInputConnection(transformPD.GetOutputPort())
		#Create a mapper and actor for the arrow
		mapper = vtk.vtkPolyDataMapper()
		actor = vtk.vtkActor()
		
		mapper.SetInputConnection(transformPD.GetOutputPort())
		
		
		actor.SetMapper(mapper)
		actor.GetProperty().SetColor(c[0],c[1],c[2])
		self.arrows.append(actor)
		self.writeObjects=True


	def addCSYS(self,orig,CSYS,scaling=1.0):
		lines=[]
		for i in range(0,3):#CSYS.shape[0]):

			ept=orig+CSYS[i]*scaling
			line = vtk.vtkLineSource()
			line.SetPoint1(orig[0],orig[1],orig[2])
			line.SetPoint2(ept[0],ept[1],ept[2])

			mapper = vtk.vtkPolyDataMapper()
			source=line
			#if vtk.VTK_MAJOR_VERSION <= 5:
				#mapper.SetInput(source.GetOutput())
			#else:
			
			#mapper.SetInputData(source)
			mapper.SetInputConnection(source.GetOutputPort())
			actor = vtk.vtkActor()
			actor.SetMapper(mapper)
			if i ==0:
				actor.GetProperty().SetColor(1,0,0)
			if i ==1:
				actor.GetProperty().SetColor(0,1,0)
			if i ==2:
				actor.GetProperty().SetColor(0,0,1)
			lines.append(actor)	
		self.axisSystems.append(lines)

	def addSphere(self,orig,radius=1.0,label='',labelOffset=np.array([1,0,0]),color=(0,0,0) ):
		source = vtk.vtkSphereSource()
		source.SetCenter(orig[0],orig[1],orig[2])
		source.SetRadius(radius)
		source.SetThetaResolution(16)
		source.SetPhiResolution(16)
		self.output.objects.AddInputConnection(source.GetOutputPort())
		mapper = vtk.vtkPolyDataMapper()
		#if vtk.VTK_MAJOR_VERSION <= 5:
			#mapper.SetInput(source.GetOutput())
		#else:
			#mapper.SetInputConnection(source.GetOutputPort())
		mapper.SetInputConnection(source.GetOutputPort())
		actor = vtk.vtkActor()
		actor.SetMapper(mapper)
		actor.GetProperty().SetColor(color[0],color[1],color[2])
		self.spheres.append(actor)
		self.addTextLabel(orig+labelOffset,label)
		

	def drawLines(self):
		for l in self.lines:
			self.ren.AddActor(l);

	def drawArrows(self):
		for l in self.arrows:
			self.ren.AddActor(l);

	def drawAxisSystems(self):
		for cs in self.axisSystems:
			for l in cs:
				self.ren.AddActor(l);

	def drawSpheres(self):
		for s in self.spheres:
			self.ren.AddActor(s);


	def addVlmPanel(self,panelI):
		if  True:#not panelI.isWake:
			cp=panelI.getCpoint()
			nv=panelI.getForce()
			f=panelI.getForce()
			ff=np.sqrt(nv[0]**2 + nv[1]**2 + nv[2]**2)
			self.output.points.InsertNextPoint(panelI.pts[0])
			self.output.points.InsertNextPoint(panelI.pts[1])
			self.output.points.InsertNextPoint(panelI.pts[2])
			self.output.points.InsertNextPoint(panelI.pts[3])

			self.output.quadCount+=1
			self.output.quadl.append(vtk.vtkQuad())
			quad=self.output.quadl[self.output.quadCount]

			qcnt=self.output.quadCount*4
			quad.GetPointIds().SetId(0,qcnt+0)
			quad.GetPointIds().SetId(1,qcnt+1)
			quad.GetPointIds().SetId(2,qcnt+2)
			quad.GetPointIds().SetId(3,qcnt+3)
			
			self.output.quads.InsertNextCell(self.output.quadl[self.output.quadCount])
	#		self.colors.InsertNextValue(ff)
			self.colors.InsertNextValue(panelI.getPressureCoef())
			self.scG.InsertNextValue(panelI.getG())
			self.scS.InsertNextValue(panelI.getS())
			self.scGid.InsertNextValue(panelI.gid)
			self.scPGid.InsertNextValue(panelI.pgid)
			
			self.force.InsertNextTuple3(f[0],f[1],f[2])
			
			if self.showAxisSystems:
				self.addCSYS(panelI.controlPoint, panelI.csys,scaling=0.01)
			self.addTextLabel(panelI.controlPoint, panelI.gid)
			self.writePolydata=True
			
	def addPolygon(self,poly,scalar=0.0):
		for i in range(0,poly.shape[0]):			
			self.output.polyVertices.InsertNextPoint(poly[i])

		polygon = vtk.vtkPolygon()
		
		polygon.GetPointIds().SetNumberOfIds(poly.shape[0]) 
		for i in range(0,poly.shape[0]):			
			polygon.GetPointIds().SetId(i, self.output.polyVerticesCount+i)

		self.output.polyVerticesCount+=poly.shape[0]

		self.output.polyList.append(polygon)
		self.output.polygons.InsertNextCell(polygon)
		self.colors.InsertNextValue(scalar)
		self.writePolydata=True
		
	def drawPolyLines(self):
		for i in range(0,len(self.pointSets)):
			cells=vtk.vtkCellArray()
			points=vtk.vtkPoints()

			ps=self.pointSets[i]
			pl=self.polylines[i]
			pl.GetPointIds().SetNumberOfIds(ps.GetNumberOfPoints())
			for j in range(0,ps.GetNumberOfPoints()):
				points.InsertNextPoint(ps.GetPoint(j))
				pl.GetPointIds().SetId(j,j)
			cells.InsertNextCell(pl)

			linePolydata = vtk.vtkPolyData()
			linePolydata.SetPoints(points)
			linePolydata.SetLines(cells)
			
			mapper = vtk.vtkPolyDataMapper();
	#		#mapper.SetLookupTable(self.lut)
	#		mapper.InterpolateScalarsBeforeMappingOn()
			#mapper.SetInput(linePolydata)
			#mapper.SetInputConnection(linePolydata.GetProducerPort())
			source=linePolydata
			if vtk.VTK_MAJOR_VERSION <= 5:
				mapper.SetInput(source.GetOutput())
			else:
				#mapper.SetInputConnection(source.GetOutputPort())
				mapper.SetInputData(source)
	#		mapper.SetScalarModeToUsePointData()
	#		mapper.ScalarVisibilityOn();
	#		mapper.SetScalarRange(self.colors.GetRange())
			surfaceActor = vtk.vtkActor();
			surfaceActor.SetMapper(mapper);
			surfaceActor.GetProperty().SetColor(1,0,0)
			self.ren.AddActor(surfaceActor);

		
	def drawQuads(self):
		#print self.colors.GetRange()
		self.output.polydata = vtk.vtkPolyData();
		self.output.polydata.SetPoints(self.output.points)
		self.output.polydata.SetPolys(self.output.quads)
		self.output.polydata.GetCellData().AddArray(self.colors)
		self.output.polydata.GetCellData().AddArray(self.scGid)
		self.output.polydata.GetCellData().AddArray(self.scPGid)
		self.output.polydata.GetCellData().AddArray(self.force)

		self.output.polydata.GetCellData().AddArray(self.scG)
		self.output.polydata.GetCellData().AddArray(self.scS)
		self.output.polydata.GetCellData().SetActiveScalars('CP')
		if vtk.VTK_MAJOR_VERSION <= 5:
			self.output.polydata.Update()	 

		self.lut=vtk.vtkLookupTable()
		self.lut.SetNumberOfTableValues(100)
		self.lut.SetTableRange(self.colors.GetRange())
		self.lut.SetHueRange(0.667, 0.0)
		self.lut.Build()

		mapper = vtk.vtkPolyDataMapper();
		mapper.SetLookupTable(self.lut)
		mapper.InterpolateScalarsBeforeMappingOn()
		#mapper.SetInput(self.output.polydata)
		#mapper.SetInputConnection(self.output.polydata.GetProducerPort())
		source=self.output.polydata
		if vtk.VTK_MAJOR_VERSION <= 5:
			mapper.SetInput(source.GetOutput())
		else:
			mapper.SetInputData(source)
#		mapper.SetScalarModeToUsePointData()
		mapper.ScalarVisibilityOn();
		mapper.SetScalarRange(self.colors.GetRange())
		surfaceActor = vtk.vtkActor();
		surfaceActor.SetMapper(mapper);
		surfaceActor.GetProperty().LightingOff()
		self.ren.AddActor(surfaceActor);


		
	def drawPolygons(self,asWire=False):
		#print self.colors.GetRange()
		self.output.polydata = vtk.vtkPolyData();
		self.output.polydata.SetPoints(self.output.polyVertices)
		self.output.polydata.SetPolys(self.output.polygons)
		self.output.polydata.GetCellData().SetScalars(self.colors)
		if vtk.VTK_MAJOR_VERSION <= 5:
			self.output.polydata.Update()	 

		self.lut=vtk.vtkLookupTable()
		self.lut.SetNumberOfTableValues(100)
		self.lut.SetTableRange(self.colors.GetRange())
		self.lut.SetHueRange(0.667, 0.0)
		self.lut.Build()

		mapper = vtk.vtkPolyDataMapper();
		mapper.SetLookupTable(self.lut)
		mapper.InterpolateScalarsBeforeMappingOn()
		#mapper.SetInput(self.output.polydata)
		#mapper.SetInputConnection(self.output.polydata.GetProducerPort())
		
		source=self.output.polydata
		if vtk.VTK_MAJOR_VERSION <=5:
			mapper.SetInput(source.GetOutput())
		else:
			mapper.SetInputData(source)
#		mapper.SetScalarModeToUsePointData()
		mapper.ScalarVisibilityOn();
		mapper.SetScalarRange(self.colors.GetRange())
		surfaceActor = vtk.vtkActor();
		surfaceActor.SetMapper(mapper);
		surfaceActor.GetProperty().LightingOff()
		if asWire:
			surfaceActor.GetProperty().SetRepresentationToWireframe()
		self.ren.AddActor(surfaceActor);


	def drawText(self):
		scl=self.textScale
		for lab in self.textLabels:
			atext = vtk.vtkVectorText()
			atext.SetText(str(lab['txt']))
			textMapper = vtk.vtkPolyDataMapper()
			textMapper.SetInputConnection(atext.GetOutputPort())
			textActor = vtk.vtkFollower()
			textActor.SetMapper(textMapper)
			textActor.SetScale(scl)
			x,y,z=lab['xyz']
			textActor.AddPosition(x,y,z)
			textActor.GetProperty().SetColor(0,0,0)

			self.ren.AddActor(textActor)
			textActor.SetCamera(self.ren.GetActiveCamera())

	def addVTPFile(self,fname):
		reader=vtk.vtkXMLPolyDataReader()
		reader.SetFileName(fname)
		reader.Update()

		self.boundBox=reader.GetOutput().GetBounds()
		self.scalarRange=reader.GetOutput().GetPointData().GetScalars().GetRange()
		self.lut=vtk.vtkLookupTable()
		self.lut.SetNumberOfTableValues(100)
		self.lut.SetTableRange(self.scalarRange)
		self.lut.SetHueRange(0.667, 0.0)
		self.lut.Build()
		mapper = vtk.vtkPolyDataMapper();
		mapper.SetLookupTable(self.lut)
#		mapper.InterpolateScalarsBeforeMappingOn()
		mapper.SetInputConnection(reader.GetOutputPort())
		mapper.SetScalarModeToUsePointData()
		mapper.ScalarVisibilityOn();
		mapper.SetScalarRange(self.scalarRange)
		surfaceActor = vtk.vtkActor();
		surfaceActor.SetMapper(mapper);

		self.ren.AddActor(surfaceActor);

	def addWidgets(self,opt=None):
		if opt!=None:
			if opt.has_key('caption'):
				self.caption=opt['caption']
			if opt.has_key('showAxes'):
				self.showAxes=opt['showAxes']
			if opt.has_key('showCompass'):
				self.showCompass=opt['showCompass']
			if opt.has_key('showScalarBar'):
				self.showScalarBar=opt['showScalarBar']
			if opt.has_key('showXYPlane'):
				self.showXYPlane=opt['showXYPlane']
			if opt.has_key('showYZPlane'):
				self.showYZPlane=opt['showYZPlane']
			if opt.has_key('showZXPlane'):
				self.showZXPlane=opt['showZXPlane']
			if opt.has_key('xlabel'):
				self.xlabel=opt['xlabel']
			if opt.has_key('ylabel'):
				self.ylabel=opt['ylabel']
			if opt.has_key('ylabel'):
				self.zlabel=opt['zlabel']
			if opt.has_key('xrange'):
				self.xrange=opt['xrange']
			if opt.has_key('yrange'):
				self.yrange=opt['yrange']
			if opt.has_key('zrange'):
				self.zrange=opt['zrange']


		prn=1000.
		pc=-prn
		plXY = vtk.vtkPlaneSource()
		plXY.SetPoint1(prn,-prn,0)
		plXY.SetPoint2(-prn,prn,0)
		plXY.SetOrigin(pc,pc,0)
		plXY.SetCenter(0,0,0)
		plXYmap = vtk.vtkPolyDataMapper()
		if vtk.VTK_MAJOR_VERSION <= 5:
			plXYmap.SetInput(plXY.GetOutput())
		else:
			plXYmap.SetInputConnection(plXY.GetOutputPort())
		
		plXYact = vtk.vtkActor()
		plXYact.SetMapper(plXYmap)
		plXYact.GetProperty().SetOpacity(0.1)

		plYZ = vtk.vtkPlaneSource()
		plYZ.SetCenter(0,pc,pc)
		plYZ.SetPoint1(0,prn,-prn)
		plYZ.SetPoint2(0,-prn,prn)
		plYZmap = vtk.vtkPolyDataMapper()
		if vtk.VTK_MAJOR_VERSION <= 5:
			plYZmap.SetInput(plYZ.GetOutput())
		else:
			plYZmap.SetInputConnection(plYZ.GetOutputPort())
		
		plYZact = vtk.vtkActor()
		plYZact.SetMapper(plYZmap)
		plYZact.GetProperty().SetOpacity(0.1)

		plZX = vtk.vtkPlaneSource()
		plZX.SetCenter(pc,0,pc)
		plZX.SetPoint1(prn,0,-prn)
		plZX.SetPoint2(-prn,0,prn)
		plZXmap = vtk.vtkPolyDataMapper()
		if vtk.VTK_MAJOR_VERSION <= 5:
			plZXmap.SetInput(plZX.GetOutput())
		else:
			plZXmap.SetInputConnection(plZX.GetOutputPort())
		
		plZXact = vtk.vtkActor()
		plZXact.SetMapper(plZXmap)
		plZXact.GetProperty().SetOpacity(0.1)

		ax=vtk.vtkAxesActor()
		ax.GetXAxisCaptionActor2D().GetProperty().SetColor(0,0,0)
		ax.GetYAxisCaptionActor2D().GetProperty().SetColor(0,0,0)
		ax.GetZAxisCaptionActor2D().GetProperty().SetColor(0,0,0)

		xa=vtk.vtkAxisActor()
		xa.SetPoint1(0,0,0)
		xa.SetPoint2(1000,0,0)
		xa.SetRange((0,1000))
		xa.SetBounds(-1.0, 1000.0, -1.0, 1.0, -1.0, 1.0)

		self.ow=vtk.vtkOrientationMarkerWidget()
		textActor = vtk.vtkTextActor()
		textActor.GetTextProperty().SetFontSize ( 22 )
		textActor.SetPosition2( 100, 100 )
		textActor.SetInput(r'3D view: '+str(self.caption)+r' Scaling: '+str(self.scalingFactor))
		textActor.GetTextProperty().SetColor ( 0.0,0.0,0.0 )
	

			
		if self.showXYPlane:
			self.ren.AddActor(plXYact)
		if self.showYZPlane:
			self.ren.AddActor(plYZact)
		if self.showZXPlane:
			self.ren.AddActor(plZXact)
		if self.showCaption:
			self.ren.AddActor2D( textActor )
		if self.showScalarBar:
			self.scalar_bar = vtk.vtkScalarBarActor()
			self.scalar_bar.SetOrientationToHorizontal()
			self.scalar_bar.SetLookupTable(self.lut)
			self.scalar_bar.SetTitle("Cf value");
			self.scalar_bar.SetNumberOfLabels(11)
			self.scalar_bar.GetProperty().SetColor ( 0.0,0.0,0.0 )
			self.scalar_bar_widget = vtk.vtkScalarBarWidget()
			self.scalar_bar_widget.SetInteractor(self.iren)
			self.scalar_bar_widget.SetScalarBarActor(self.scalar_bar)
			self.scalar_bar_widget.On()
		if self.showCompass:
			self.ow.SetOrientationMarker(ax)
			self.ow.SetInteractor(self.iren)
#			self.ow.SetViewport( 0.0, 0.0, 0.4, 0.4 )
			self.ow.SetEnabled( 1 )
			self.ow.InteractiveOn()
		if self.showAxes:
			c=vtk.vtkCubeAxesActor()
			c.SetBounds(self.boundBox)
			c.SetXTitle(self.xlabel)
			c.SetYTitle(self.ylabel)
			c.SetZTitle(self.zlabel)

			if self.xrange != None:
				c.SetXAxisRange(self.xrange)
			if self.yrange != None:
				c.SetYAxisRange(self.yrange)
			if self.zrange != None:
				c.SetZAxisRange(self.zrange)
			c.GetProperty().SetColor(0., 0., 0.)
			c.SetCamera(self.ren.GetActiveCamera())
			self.ren.AddActor(c)

