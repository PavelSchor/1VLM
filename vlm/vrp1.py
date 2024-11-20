from panelMethod.vlmPanels import VlmPanel
from panelMethod.vortexPanels import VortexRingPanel
import numpy as np
import copy
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
		self.nLPanels=-1
		self.LPanels={}
		self.refPT=np.zeros(3)

	def addPanel(self,panel):
		self.nPanels=self.nPanels+1
		self.panels[self.nPanels]=panel
		if not panel.isWake:
			self.nLPanels+=1;
			self.LPanels[self.nLPanels]=panel

	def setRefPT(self,PT):
		self.refPT=PT

	def getPanelList(self):
		return self.panels.keys()

	def getNumberOfPanels(self):
		return len(self.panels)

	def getPanel(self,i):
		return self.panels[i]

	def getLiftPanelList(self):
		return self.LPanels.keys()

	def getNumberOfLiftPanels(self):
		return len(self.LPanels)

	def getLiftPanel(self,i):
		return self.LPanels[i]


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

	def addSubsection(self,r):
		for i in r.getPanelList():
			self.addPanel(r.getPanel(i))

	def addPanelGroup(self,pg):
		self.addSubsection(pg)

class GlobalDomain(GlobalEnvironment):
	def __init__(self):
		self.refPT=np.zeros(3)
		self.regions={}
		self.gids={}
		self.nReg=-1
		self.nPan=0
		self.nLPan=0
		self.panelCounter=0
		
	def panelCounterIncrease(self):
		self.panelCounter+=1

	def buildGids(self):
		for r in self.getRegionList():
			for i in self.getRegion(r).getPanelList():
				
				self.getRegion(r).getPanel(i).connectFlowDomain(self)
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
		self.nLPan+=self.regions[self.nReg].getNumberOfLiftPanels()
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
		vel=copy.copy(self.getFreeVelocity())
		for r in self.getRegionList():
			for i in self.getRegion(r).getPanelList():
				vel+=self.getRegion(r).getPanel(i).getUVW(PT)
		return vel

	def getPHIAtPoint(self,PT):
		phi=0.0
		for r in self.getRegionList():
			for i in self.getRegion(r).getPanelList():
				phi+=self.getRegion(r).getPanel(i).getPHI(PT)
		return phi

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

	def getPanelOffset(self,i):
		ofs=0
		for i in range(0,i):
			ofs=ofs+self.getRegion(i).getNumberOfPanels()
		return ofs		

	def getLiftPanelOffset(self,i):
		ofs=0
		for i in range(0,i):
			ofs=ofs+self.getRegion(i).getNumberOfLiftPanels()
		return ofs		


	def solveAij(self):
		self.buildGids()
		print 'assembling Aij matrix'
		print strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
		vel=self.getFreeVelocity()
		A=np.zeros((self.nLPan,self.nLPan))
		B=np.zeros((self.nLPan,1))
		b=np.zeros((self.nLPan,self.nLPan))
		s=np.zeros((self.nLPan,1))		

		for r1 in self.getRegionList():
			for r2 in self.getRegionList():
				ofsI=self.getLiftPanelOffset(r1)
				ofsJ=self.getLiftPanelOffset(r2)
				nPanI=self.getRegion(r1).getNumberOfLiftPanels()
				nPanJ=self.getRegion(r2).getNumberOfLiftPanels()
				for i in range(0,nPanI):
					panelI=self.getRegion(r1).getLiftPanel(i)
					panelI.setS(1.0)#np.dot(self.getFreeVelocity(),panelI.getNV()))
					s[i]=np.dot(self.getFreeVelocity(),panelI.normalVect)
					for j in range(0,nPanJ):
						panelJ=self.getRegion(r2).getLiftPanel(j)
						b[i+ofsI][j+ofsJ]=np.dot( panelI.normalVect , panelJ.getSourceUVW(panelI.getCpoint(),1.0))

		self.b=b
		self.s=s
		for r1 in self.getRegionList():
			for r2 in self.getRegionList():
				ofsI=self.getLiftPanelOffset(r1)
				ofsJ=self.getLiftPanelOffset(r2)
				nPanI=self.getRegion(r1).getNumberOfLiftPanels()
				nPanJ=self.getRegion(r2).getNumberOfLiftPanels()
				for i in range(0,nPanI):
					panelI=self.getRegion(r1).getLiftPanel(i)
					for j in range(0,nPanJ):
						panelJ=self.getRegion(r2).getLiftPanel(j)
						vij=panelJ.getUVW(panelI.getCpoint())
						if not  panelI.isWake:
							A[i+ofsI][j+ofsJ]=np.dot(vij,panelI.getNV()) 
						
						wakeInf=0.0
						if  panelJ.wakes:
							for wp in range(0,len(panelJ.wakes)):
								wPan=self.getPanelByGid(panelJ.wakes[wp])
								vwake=wPan.getUVW(panelI.getCpoint())
								wakeInf+=np.dot(vwake,panelI.getNV() )
#						wakeInf=0.0
						if panelJ.isLast:
							A[i+ofsI][j+ofsJ]+=wakeInf
						if panelJ.isFirst:
							A[i+ofsI][j+ofsJ]-=wakeInf
						B[i+ofsI]-=np.dot(panelJ.getSourceUVW(panelI.getCpoint(),panelJ.S) ,panelI.getNV())
#				B[i+ofsI]-=np.dot(vel,panelI.getNV())
	
				
		self.B=B#(1.)*np.dot(b,s)
		self.B0=B
		self.A=A
		print 'Aij matrix assembled'
		print 'solving Aij'
		G=np.squeeze(np.linalg.solve(A,self.B))
#		G=np.squeeze(np.linalg.lstsq(A,self.B)[0])
		print 'postprocessing results'
		print strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
		self.G=G
#		self.B0=B0
		for r in self.getRegionList():
			nPan=self.getRegion(r).getNumberOfLiftPanels()
			ofs=self.getLiftPanelOffset(r)	
			for i in range(0,nPan):
				self.getRegion(r).getLiftPanel(i).setG(G[ofs+i])
				#########################
				p=self.getRegion(r).getLiftPanel(i)
				if p.isLast:
					if p.wakes:
						g= self.getPanelByGid(p.firstPanel).G - p.G
						for ppw in range(0,len(p.wakes)):
							self.getPanelByGid(p.wakes[ppw]).setG(g)
#							print "panel: "+str(self.getPanelByGid(p.wakes[ppw]).gid)+" g is "+str(g)
		print 'panels G set'
		print strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())			
		for r in self.getRegionList():
			nPan=self.getRegion(r).getNumberOfLiftPanels()
			ofs=self.getLiftPanelOffset(r)
			for i in range(0,nPan):	
				panelI=self.getRegion(r).getLiftPanel(i)
				CPT=panelI.getCpoint()
				velCPT=self.getVelocityAtPoint(CPT)
				panelI.setFreeStream(velCPT)
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
		dom1.setPanelsG(0.0)
		dom1.solveAij()

	def plotResults(self,lineScaling=1.0):
		self.view=PanelViewer3D()
		dom1=self.dom1
		for r in dom1.getRegionList():
			nPan=dom1.getRegion(r).getNumberOfPanels()
			lc=self.view.newPolyLine()
			for i in range(0,nPan):
				panelI=dom1.getRegion(r).getPanel(i)
				if not panelI.isWake:
					npt=panelI.getCpoint()+ np.array([0,1,0])* panelI.getPressureCoef()*lineScaling
					#panelI.getCpoint()+panelI.getNV()*panelI.getPressureCoef()*lineScaling
					if (i==0):
						self.view.pointSets[lc].InsertNextPoint(panelI.getCpoint())
					self.view.pointSets[lc].InsertNextPoint(npt)
					if (i==nPan-1):
						self.view.pointSets[lc].InsertNextPoint(panelI.getCpoint())
					self.view.addVlmPanel(panelI)
					
		self.view.drawQuads()
		self.view.drawPolyLines()
		self.view.addWidgets({'showAxes':False,'showScalarBar':True,'showXYPlane':False,'showYZPlane':False,'showZXPlane':False})
		self.view.show()
