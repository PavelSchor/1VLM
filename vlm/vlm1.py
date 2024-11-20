from panelMethod.vlmPanels import VlmPanel
from panelMethod.vortexPanels import VortexRingPanel
import numpy as np
from panelMethod.env import GlobalEnvironment
from utils.transform import EulerMatrix
from utils.b3Vect import B3vect 
vect=B3vect()
from panelMethod.viewer import *

from time import gmtime, strftime


class subSection(GlobalEnvironment):
	def __init__(self):
		self.nPanels=-1
		self.panels={}
		self.refPT=np.zeros(3)

	def addPanel(self,panel):
		self.nPanels=self.nPanels+1
		self.panels[self.nPanels]=panel

	def setRefPT(self,PT):
		self.refPT=PT

	def getPanelList(self):
		return self.panels.keys()

	def getNumberOfPanels(self):
		return len(self.panels)


	def getPanel(self,i):
		return self.panels[i]

	def getForceAtPT(self,PT):
		f=np.zeros(3)
		m=np.zeros(3)
		for i in self.getPanelList():
			panelI=self.getPanel(i)
			f=f+panelI.getForce()
			r=panelI.getVortexPoint0()-PT
			m=m+np.cross(panelI.getForce(),r)
		return np.hstack((f,m))

	def getForceAtRefPT(self):
		return self.getForceAtPT(self.refPT)

	def getArea(self):
		A=0.0
		for i in self.getPanelList():
			A+=self.getPanel(i).getArea()
		return A
    
class PanelGroup(subSection):
	def __init__(self):
		self.nPanels=-1
		self.panels={}
		self.refPT=np.zeros(3)

	def addPanel(self,panel):
		self.nPanels=self.nPanels+1
		self.panels[self.nPanels]=panel
#		print panel
#		print self.panels[self.nPanels]

	def addSubsection(self,r):
		for i in r.getPanelList():
			self.addPanel(r.getPanel(i))

	def addPanelGroup(self,pg):
		self.addSubsection(pg)
    
    def getSubsectionByMask(self,msk):
        pg=PanelGroup()
        for i in range(0,len(self.panels)):
            if msk[i]:
                pg.addPanel(self.panels[i])
        return pg
        
        
class GlobalDomain(GlobalEnvironment):
	def __init__(self):
		self.refPT=np.zeros(3)
		self.regions={}
		self.gids={}
		self.nReg=-1
		self.nPan=0
		self.panelCounter=0
		
	def panelCounterIncrease(self):
		self.panelCounter+=1

	def buildGids(self):
		for r in self.getRegionList():
			for i in self.getRegion(r).getPanelList():
				gid=self.getRegion(r).getPanel(i).gid
				self.gids[gid]=self.getRegion(r).getPanel(i)

	def getPanelByGid(self,gid):
		return self.gids[gid]
				
	def setRefPT(self,PT):
		self.refPT=PT

	def addRegion(self,reg):
		self.nReg=self.nReg+1
		self.regions[self.nReg]=reg
		self.nPan=self.nPan+self.regions[self.nReg].getNumberOfPanels()

#		print '===================== add region ================================='
#		print reg 
#		print self.regions[self.nReg]
#
#		for i in range(0,len(reg.panels)):
#			print reg.panels[i]
#			print self.regions[self.nReg].panels[i]
	
	def getRefPT(self):
		return self.refPT

	def getRegionList(self):
		return self.regions.keys()

	def getRegion(self,i):
		return self.regions[i]

	def getVelocityAtPoint(self,PT):
		vel=self.getFreeVelocity()
		for r in self.getRegionList():
		#	nPan=self.getRegion(r).getNumberOfPanels()
		#	ofs=self.getPanelOffset(r)	
		#	for i in range(0,nPan):
			for i in self.getRegion(r).getPanelList():
				panelI=self.getRegion(r).getPanel(i)	
				vel=vel+panelI.getUVW(PT)
		return vel

	def reinitPanelsGeometry(self):
		for r in self.getRegionList():
			for i in self.getRegion(r).getPanelList():
				self.getRegion(r).getPanel(i).initGeometry()

	def setPanelsG(self,G):
		for r in self.getRegionList():
			for i in self.getRegion(r).getPanelList():
				self.getRegion(r).getPanel(i).setG(G)

	def setPanelsFreeStream(self,v):
		for r in self.getRegionList():
			for i in self.getRegion(r).getPanelList():
				self.getRegion(r).getPanel(i).setFreeStream(v)
			
	def getDownwashAtPoint(self,PT):
		vel=np.zeros(3)
		for r in self.getRegionList():
			nPan=self.getRegion(r).getNumberOfPanels()
			for i in range(0,nPan):
				vel=vel+self.getRegion(r).getPanel(i).getUVW(PT,sf0=0.0,sf1=1.0,sf2=1.0)
		return vel


	def getPanelOffset(self,i):
		ofs=0
		for i in range(0,i):
			ofs=ofs+self.getRegion(i).getNumberOfPanels()
		return ofs		

	def solveAij(self):
		print 'assembling Aij matrix'
		print strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
		vel=self.getFreeVelocity()
		A=np.zeros((self.nPan,self.nPan))
		B=np.zeros((self.nPan,1))
		for r1 in self.getRegionList():
			for r2 in self.getRegionList():
				ofsI=self.getPanelOffset(r1)
				ofsJ=self.getPanelOffset(r2)
				nPanI=self.getRegion(r1).getNumberOfPanels()
				nPanJ=self.getRegion(r2).getNumberOfPanels()
				for i in range(0,nPanI):
					panelI=self.getRegion(r1).getPanel(i)
					for j in range(0,nPanJ):
						panelJ=self.getRegion(r2).getPanel(j)
						vij=panelJ.getUVW(panelI.getCpoint())
						A[i+ofsI][j+ofsJ]=np.dot(vij,panelI.getNV()) 
					B[i+ofsI]=-np.dot(vel,panelI.getNV())

		print 'Aij matrix assembled'
		print 'solving Aij'
		S=np.squeeze(np.linalg.solve(A,B))
#		S=np.squeeze(np.linalg.lstsq(A,B)[0])
		print 'postprocessing results'
		print strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
		self.A=A
		self.S=S
		self.B=B
		for r in self.getRegionList():
			nPan=self.getRegion(r).getNumberOfPanels()
			ofs=self.getPanelOffset(r)	
			for i in range(0,nPan):
				panelI=self.getRegion(r).getPanel(i)	
				panelI.setG(S[ofs+i])
		print 'panels G set'
		print strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())			
		for r in self.getRegionList():
			nPan=self.getRegion(r).getNumberOfPanels()
			ofs=self.getPanelOffset(r)
			for i in range(0,nPan):	
				panelI=self.getRegion(r).getPanel(i)
				CPT=panelI.getCpoint()
				velCPT=self.getVelocityAtPoint(CPT)
				panelI.setFreeStream(velCPT)
				#panelI.dw=np.zeros(3)#self.getDownwashAtPoint(CPT)
				panelI.setForce()
		print 'results computed'
		print strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
class VlmProblem(GlobalEnvironment):
	def __init__(self):
		self.vel=np.zeros(3)
		self.velMagn=10.0
		self.alpha=0.0
		self.beta=0.0
		self.pressure=1.0	
		self.dom1=GlobalDomain()

	def setAlpha(self,a):
		self.alpha=a
	
	def setBeta(self,b):
		self.beta=b

	def setVelocityMagn(self,v):
		self.velMagn=v
		vel=np.array([v,0.,0.])
		a=self.getAlpha()
		b=self.getBeta()
		eMatrix=EulerMatrix(0.0,a,b)
		self.setVelocity(eMatrix.transform(vel))
	
	def getAlpha(self):
		return self.alpha

	def getBeta(self):
		return self.beta

	def setVelocity(self,vel):
		GlobalEnvironment.setFreeVelocity(vel)

	def getVelocity(self):
		return self.getFreeVelocity()
	
	def addRegion(self,reg):
		self.dom1.addRegion(reg)
		
	def getSolution(self):
		dom1=self.dom1
		dom1.setPanelsFreeStream(self.getFreeVelocity())
		dom1.reinitPanelsGeometry()
		dom1.setPanelsG(1.0)
		dom1.solveAij()

	def plotResults(self,lineScaling=1.0):
		self.view=PanelViewer3D()
		dom1=self.dom1
		for r in dom1.getRegionList():
			nPan=dom1.getRegion(r).getNumberOfPanels()
			lc=self.view.newPolyLine()
			for i in range(0,nPan):
				panelI=dom1.getRegion(r).getPanel(i)
				np=panelI.getVortexPoint0()+panelI.getNV()*panelI.getPressureCoef()*lineScaling
				if (i==0):
					self.view.pointSets[lc].InsertNextPoint(panelI.getVortexPoint0())
				self.view.pointSets[lc].InsertNextPoint(np)
				if (i==nPan-1):
					self.view.pointSets[lc].InsertNextPoint(panelI.getVortexPoint0())
				self.view.addVlmPanel(panelI)
				
		self.view.drawQuads()
		self.view.drawPolyLines()
		self.view.addWidgets({'showAxes':False,'showScalarBar':True,'showXYPlane':False,'showYZPlane':False,'showZXPlane':False})
		self.view.show()
