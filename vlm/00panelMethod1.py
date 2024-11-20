from panelMethod.doubletSourcePanels import doubletSourcePanel
import numpy as np
import copy
from panelMethod.env import GlobalEnvironment
from utils.transform import EulerMatrix
from utils.b3Vect import B3vect 
vect=B3vect()
from panelMethod.viewer import *

from time import gmtime, strftime
from scipy.sparse.linalg import gmres

import ctypes
from numpy.ctypeslib import ndpointer

from Cython import nogil

try:
	testlib = ctypes.CDLL('panelMethod/testlib_w32.so')
except:
	pass
try:
	testlib = ctypes.CDLL('panelMethod/testlib_w64.so')
except:
	pass


try:
	testlib = ctypes.CDLL('panelMethod/testlib_lnx32.so')
except:
	pass

try:
	testlib = ctypes.CDLL('panelMethod/testlib_lnx64.so')
except:
	pass

#void getPHI_at_points(double*PHI, double*poi, int n, double*pts,int npt, int ptdim_max,int*ptdim, double*ORIGS, double*csysesLoc,double*D,double*S){
testlib.getPHI_at_points.restype = ctypes.c_void_p;
testlib.getPHI_at_points.argtypes = [ndpointer(ctypes.c_double),
ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
ctypes.c_int,
ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
ctypes.c_int,
ctypes.c_int,
ndpointer(ctypes.c_int),
ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
ndpointer(ctypes.c_double),
ndpointer(ctypes.c_double)];

#void getAij_matrix(double*AIJ, double*poi, int n, double*pts,int npt, int ptdim_max,int*ptdim, double*ORIGS, double*csysesLoc,double*D,double*S,int*gid,int*isFirst,int*isLast,int*isWake){

testlib.getAij_matrix.restype = ctypes.c_void_p;
testlib.getAij_matrix.argtypes = [ndpointer(ctypes.c_double),
ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
ctypes.c_int,
ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
ctypes.c_int,
ctypes.c_int,
ndpointer(ctypes.c_int),
ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
ndpointer(ctypes.c_double),
ndpointer(ctypes.c_double),
ndpointer(ctypes.c_int),
ndpointer(ctypes.c_int),
ndpointer(ctypes.c_int),
ndpointer(ctypes.c_int),
ndpointer(ctypes.c_int),
ndpointer(ctypes.c_int)];




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
			r=panelI.midPoint -PT
			m=m+np.cross(r,panelI.getForce())
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
		self.view=PanelViewer3D()
		self.refPT=np.zeros(3)
		self.regions={}
		self.gids={}
		self.nReg=-1
		self.nPan=0
		self.nLPan=0
		self.panelCounter=0
		self.wakePanels=[]
		
	def panelCounterIncrease(self):
		self.panelCounter+=1

	def buildGids(self):
		for r in self.getRegionList():
			for i in self.getRegion(r).getPanelList():
				
				self.getRegion(r).getPanel(i).connectFlowDomain(self)
				gid=self.getRegion(r).getPanel(i).gid				
				self.gids[gid]=self.getRegion(r).getPanel(i)
				if self.getRegion(r).getPanel(i).isWake:
					self.wakePanels.append(gid) 
		
	def getPanelByGid(self,gid):
		return self.gids[gid]
				
	def setRefPT(self,PT):
		self.refPT=PT

	def addRegion(self,reg):
		self.nReg=self.nReg+1
		self.regions[self.nReg]=reg
		self.nPan=self.nPan+self.regions[self.nReg].getNumberOfPanels()
		self.nLPan+=self.regions[self.nReg].getNumberOfLiftPanels()
		self.buildGids()
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
		vel=np.zeros(3)
		vel+=self.getFreeVelocity()
		for r in self.getRegionList():
			for i in self.getRegion(r).getPanelList():
				vel+=self.getRegion(r).getPanel(i).getUVW(PT)
		return vel

	def getPHI0AtPoint(self,PT):
		return np.dot(self.getFreeVelocity(),PT)



	def getPHIAtPoint(self,PT):
		phi=self.getPHI0AtPoint(PT)
	#	phi=0.0
		for r in self.getRegionList():
			for i in self.getRegion(r).getPanelList():
				phi=phi+self.getRegion(r).getPanel(i).getPHI(PT)
		return phi
	def getPHIInternalAtPoint(self,PT):
		phi=0.0
		for r in self.getRegionList():
			for i in self.getRegion(r).getPanelList():
				phi=phi+self.getRegion(r).getPanel(i).getPHI(PT)
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


	def solveAij(self,normalsOrientation=1.0):
		self.buildGids()
		print 'assembling Aij matrix'
		print strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())

####
####		for r in self.getRegionList():
####			nPan=self.getRegion(r).getNumberOfLiftPanels()
####			ofs=self.getLiftPanelOffset(r)	
####			for i in range(0,nPan):
####				p=self.getRegion(r).getLiftPanel(i)
#####				self.getRegion(r).getLiftPanel(i).setS(normalsOrientation*np.dot(p.normalVect,self.getFreeVelocity()))
####				vloc=p.transform2loc(np.zeros(3),self.getFreeVelocity())
####				self.getRegion(r).getLiftPanel(i).setS(normalsOrientation * vloc[2])
####		for r1 in self.getRegionList():
####			for r2 in self.getRegionList():
####				ofsI=self.getPanelOffset(r1)
####				nPanI=self.getRegion(r1).getNumberOfPanels()
####				for i in range(0,nPanI):
####					panelI=self.getRegion(r1).getPanel(i)
####					panelI.G=0.0
####
####
####		A=np.zeros((self.nPan,self.nPan))
####		B=np.zeros(self.nPan)
####
####		for r1 in self.getRegionList():
####			for r2 in self.getRegionList():
####				ofsI=self.getPanelOffset(r1)
####				nPanI=self.getRegion(r1).getNumberOfPanels()
####				for i in range(0,nPanI):
####					panelI=self.getRegion(r1).getPanel(i)
####					B[i+ofsI]=self.getPHIInternalAtPoint(panelI.controlPoint)
####
####
####		for r1 in self.getRegionList():
####			for r2 in self.getRegionList():
####				ofsI=self.getPanelOffset(r1)
####				ofsJ=self.getPanelOffset(r2)
####				nPanI=self.getRegion(r1).getNumberOfPanels()
####				nPanJ=self.getRegion(r2).getNumberOfPanels()
####				for i in range(0,nPanI):
####					panelI=self.getRegion(r1).getPanel(i)
####					for j in range(0,nPanJ):
####						panelJ=self.getRegion(r2).getPanel(j)
####						if not panelI.isWake:
####							A[i+ofsI][j+ofsJ]=panelJ.getDoubletPHI(panelI.getCpoint(),1.0) 
####						else:
####							if panelJ.isWake:
####								A[j+ofsJ][j+ofsJ]=1.0
####								
####							if panelJ.isFirst:
####								if panelI.firstPanel == panelJ.gid:
####									A[i+ofsI][j+ofsJ]=1.0
####
####							if panelJ.isLast:
####								if panelI.lastPanel == panelJ.gid:
####									A[i+ofsI][j+ofsJ]=-1.0
####
####					if panelI.isWake:
####						B[i+ofsI]=0.0
####
#######################################################################################################################################				
####		self.B=-B#(1.)*np.dot(self.b,self.d)
#####		self.B= -np.dot(self.s,self.b)
####		self.B0=B
####		self.A=A
####		print 'Aij matrix assembled'
####		print 'solving Aij'
####
####		self.G=gmres(self.A,self.B,tol=1e-16)[0]
#####		self.G=np.squeeze(np.linalg.solve(A,self.B))
#####		self.G=np.squeeze(np.linalg.lstsq(A,self.B)[0])
####		print 'postprocessing results'
####		print strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
####
#####		for r in self.getRegionList():
#####			nPan=self.getRegion(r).getNumberOfLiftPanels()
#####			ofs=self.getLiftPanelOffset(r)	
#####			for i in range(0,nPan):
#####				self.getRegion(r).getPanel(i).setG(self.G[ofs+i])
#####				p=self.getRegion(r).getLiftPanel(i)
#####				if p.isLast:
#####					if p.wakes:
#####						g= p.G - self.getPanelByGid(p.firstPanel).G #- p.G
#####						for ppw in range(0,len(p.wakes)):
#####							self.getPanelByGid(p.wakes[ppw]).setG(g)
#####
####		for r in self.getRegionList():
####			nPan=self.getRegion(r).getNumberOfPanels()
####			ofs=self.getPanelOffset(r)	
####			for i in range(0,nPan):
####				self.getRegion(r).getPanel(i).setG(self.G[ofs+i])
####				#########################
####
####
####
####		print 'panels G set'
####		print strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())			
####		for r in self.getRegionList():
####			nPan=self.getRegion(r).getNumberOfPanels()
####			ofs=self.getPanelOffset(r)
####			for i in range(0,nPan):	
####				panelI=self.getRegion(r).getPanel(i)
####				CPT=panelI.getCpoint()
####		#		panelI.phiInt=self.getPHIInternalAtPoint(CPT)
#####				panelI.setFreeStream(self.getVelocityAtPoint(CPT))
####		for r in self.getRegionList():
####			nPan=self.getRegion(r).getNumberOfPanels()
####			ofs=self.getPanelOffset(r)
####			for i in range(0,nPan):	
####				panelI=self.getRegion(r).getPanel(i)
####				panelI.setForce()
####		print 'results computed'
####		print strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
###### 
#void getPHI_at_points(double*PHI, double*poi, int n, double*pts,int npt, int ptdim_max,int*ptdim, double*ORIGS, double*csysesLoc,double*D,double*S){
#		testlib.getPHI_at_points.restype = ctypes.c_void_p;
#		testlib.getPHI_at_points.argtypes = [ndpointer(ctypes.c_double),
#		ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
#		ctypes.c_int,
#		ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
#		ctypes.c_int,
#		ctypes.c_int,
#		ndpointer(ctypes.c_int),
#		ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
#		ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
#		ndpointer(ctypes.c_double),
#		ndpointer(ctypes.c_double)];	

		self.B=np.zeros(self.nPan,dtype=np.float64)
		poi=np.zeros((self.nPan,3),dtype=np.float64)
		n=self.nPan
		pts=np.zeros((self.nPan,3*4),dtype=np.float64)
		npt=self.nPan
		ptdim_max=4
		ptdim=np.ones(self.nPan,dtype=np.int32)*4
		origs=np.zeros((self.nPan,3),dtype=np.float64)
		csysesLoc=np.zeros((self.nPan,9),dtype=np.float64)
		D=np.zeros(self.nPan)
		S=np.zeros(self.nPan)
		print 'assembling RHS by testlib '+strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
		for r in self.getRegionList():
			nPan=self.getRegion(r).getNumberOfPanels()
			ofsI=self.getPanelOffset(r)
			for i in range(0,nPan):	
				p=self.getRegion(r).getPanel(i)
				poi[i+ofsI]=np.float64(p.controlPoint)
				pts[i+ofsI]=np.float64( p.pts.flatten())
				origs[i+ofsI]=np.float64(p.controlPoint)
				csysesLoc[i+ofsI]=np.float64(p.csys.flatten())
				D[i+ofsI]=np.float64(0.0)
				vloc=p.transform2loc(np.zeros(3),self.getFreeVelocity())
				if p.isWake:
					S[i+ofsI]=np.float64(0.0)
				else:
					S[i+ofsI]=np.float64(normalsOrientation * vloc[2])
		with nogil:
			testlib.getPHI_at_points(self.B,poi,n,pts,npt,ptdim_max,ptdim,origs,csysesLoc,D,S)
		for r in self.getRegionList():
			nPan=self.getRegion(r).getNumberOfPanels()
			ofsI=self.getPanelOffset(r)
			for i in range(0,nPan):	
				p=self.getRegion(r).getPanel(i)
				if p.isWake:
					self.B[i+ofsI]=0.0	
#		self.poi, self.pts, self.origs, self.csysesLoc, self.S = poi, pts,origs, csysesLoc, S
		del poi, pts,origs, csysesLoc,D,S
#void getAij_matrix(double*AIJ, double*poi, int n, double*pts,int npt, int ptdim_max,int*ptdim, double*ORIGS, double*csysesLoc,double*D,double*S,int*gid,int*isFirst,int*isLast,int*isWake)

		self.A=np.zeros(self.nPan*self.nPan ,dtype=np.float64)
		poi=np.zeros((self.nPan,3),dtype=np.float64)
		n=self.nPan#*self.nPan
		pts=np.zeros((self.nPan,3*4),dtype=np.float64)
		npt=self.nPan
		ptdim_max=4
		ptdim=np.ones(self.nPan,dtype=np.int32)*4
		origs=np.zeros((self.nPan,3),dtype=np.float64)
		csysesLoc=np.zeros((self.nPan,9),dtype=np.float64)
		D=np.zeros(self.nPan,dtype=np.float64)
		S=np.zeros(self.nPan,dtype=np.float64)
		gid=np.zeros(self.nPan,dtype=np.int32)
		isFirst=np.ones(self.nPan,dtype=np.int32)*-1
		isLast= np.ones(self.nPan,dtype=np.int32)*-1
		isWake= np.ones(self.nPan,dtype=np.int32)*-1
		firstPanel=np.ones(self.nPan,dtype=np.int32)*-1
		lastPanel=np.ones(self.nPan,dtype=np.int32)*-1
		print 'assembling AIJ by testlib '+strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
		for r in self.getRegionList():
			nPan=self.getRegion(r).getNumberOfPanels()
			ofsI=self.getPanelOffset(r)
			for i in range(0,nPan):	
				p=self.getRegion(r).getPanel(i)

				poi[i+ofsI]=np.float64(p.controlPoint)
				pts[i+ofsI]=np.float64( p.pts.flatten())
				origs[i+ofsI]=np.float64(p.controlPoint)
				csysesLoc[i+ofsI]=np.float64(p.csys.flatten())
				D[i+ofsI]=np.float64(1.0)
				S[i+ofsI]=np.float64(0.0)

				gid[i+ofsI]=np.int32(p.gid)
				if p.isFirst:
					isFirst[i+ofsI]=np.int32(1)
				if p.isLast:
					isLast[i+ofsI]=np.int32(1)
				if p.isWake:
					isWake[i+ofsI]=np.int32(1)
				if p.firstPanel != None:
					firstPanel[i+ofsI]=np.int32(p.firstPanel)
				if p.lastPanel != None:
					lastPanel[i+ofsI]=np.int32(p.lastPanel)

		testlib.getAij_matrix(self.A,poi,n,pts,npt,ptdim_max,ptdim,origs,csysesLoc,D,S,gid,isFirst,isLast,isWake,firstPanel,lastPanel)
		del poi,pts,ptdim,origs,csysesLoc, D, S, gid,isFirst, isLast, isWake, firstPanel, lastPanel
		
		self.A=self.A.reshape((self.nPan,self.nPan))
#		self.islast=isLast
#		self.lastPanel=lastPanel
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
						if  panelI.isWake:
							if panelJ.isWake:
								self.A[j+ofsJ][j+ofsJ]=1.0
								
							if panelJ.isFirst:
								if panelI.firstPanel == panelJ.gid:
									self.A[i+ofsI][j+ofsJ]=1.0

							if panelJ.isLast:
								if panelI.lastPanel == panelJ.gid:
									self.A[i+ofsI][j+ofsJ]=-1.0

		self.B=-self.B#(1.)*np.dot(self.b,self.d)
		print 'Aij matrix assembled'
		print 'solving Aij'

		self.G=gmres(self.A,self.B,tol=1e-16)[0]
		print 'postprocessing results'
		print strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())

		for r in self.getRegionList():
			nPan=self.getRegion(r).getNumberOfPanels()
			ofs=self.getPanelOffset(r)	
			for i in range(0,nPan):
				self.getRegion(r).getPanel(i).setG(self.G[ofs+i])

		print 'panels G set'
		print strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())			
		for r in self.getRegionList():
			nPan=self.getRegion(r).getNumberOfPanels()
			ofs=self.getPanelOffset(r)
			for i in range(0,nPan):	
				panelI=self.getRegion(r).getPanel(i)
				CPT=panelI.getCpoint()
		for r in self.getRegionList():
			nPan=self.getRegion(r).getNumberOfPanels()
			ofs=self.getPanelOffset(r)
			for i in range(0,nPan):	
				panelI=self.getRegion(r).getPanel(i)
				panelI.setForce()
		print 'results computed'
		print strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())



class VlmProblem(GlobalEnvironment):
	def __init__(self):
		self.vel=np.zeros(3)
		self.velMagn=1.0
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

	def getSolutionAeroelastic(self,stopF,stopDim=None,fullRun=False,err=1e-6,maxIt=10):
		if stopDim==None:
			stop=stopF
		ii=np.linspace(0,maxIt-1,maxIt)
		fi=np.zeros(maxIt)
		for i in range(0,maxIt):
			print('solving aeroelastic iteration '+str(i))
			if stopDim==None:
				stop=stopF()
			else:
				stop=stopF()[stopDim]
			stopOld=copy.copy(stop)
			force0=copy.copy(self.elastic.force2Fem())
#			print force0
			self.elastic.resetGeometry()
			self.fem.force=force0
			#fem.constrains=constrains
			self.fem.processInputs()
			self.fem.solveDisplacements()
			self.elastic.deform(self.fem.NODE,self.fem.D)
			self.elastic.buildCsys()
			self.getSolution()
			if stopDim==None:
				stop=stopF()
			else:
				stop=stopF()[stopDim]
			fi[i]=stop
			print 'stop: '+str(stop)+' stopOld: '+str(stopOld)
			if np.abs(stop - stopOld) < err and not fullRun:
				break
			self.elastic.buildCsys()
		return ii,fi
#			print('Force: '+str(wing.getForceAtRefPT()))
#	print fem.F.toarray()[6:42]
#			f0=wing.getForceAtRefPT()
#			f=f0#/0.5/rho/sp**2.0/area
#			cm=f[5]
#			cl=f[1]
#			cd=f[0]
#			cl=-f[0]*np.sin(alpha)+f[1]*np.cos(alpha)
#			cd= f[0]*np.cos(alpha)+f[1]*np.sin(alpha)
#			cm=f[5]
#			fi[i]=fem.D[0][1]


	def plotResults(self,redraw=True,lineScaling=1.0):
		if redraw:
			self.view=PanelViewer3D()
		dom1=self.dom1
		for r in dom1.getRegionList():
			nPan=dom1.getRegion(r).getNumberOfPanels()
			lc=self.view.newPolyLine()
			for i in range(0,nPan):
				panelI=dom1.getRegion(r).getPanel(i)
				if not panelI.isWake:
					pp=panelI.getCpoint()
					pxz=np.array([pp[0],0.0,pp[2]])
#					npt=pxz+ np.array([0,-1,0])* panelI.getPressureCoef()*lineScaling
					npt=panelI.getCpoint()+ np.array([0,1,0])* panelI.getPressureCoef()*lineScaling
#					npt=panelI.getCpoint()+panelI.getNV()*panelI.getPressureCoef()*lineScaling
					nv=panelI.getCpoint()+panelI.forceVector*lineScaling
					self.view.addLine(pp,nv)
					if (i==0):
						self.view.pointSets[lc].InsertNextPoint(panelI.getCpoint())
					self.view.pointSets[lc].InsertNextPoint(npt)
					if (i==nPan-1):
						self.view.pointSets[lc].InsertNextPoint(panelI.getCpoint())
					self.view.addVlmPanel(panelI)
		try:
			for i in self.elastic.interactors:
				I=self.elastic.interactors[i]
				self.view.addCSYS(I.orig,I.csysLoc,scaling=0.2)
		except:
			pass

		self.view.drawLines()
		self.view.drawAxisSystems()
		self.view.drawQuads()
#		self.view.drawPolyLines()
		self.view.addWidgets({'showAxes':False,'showScalarBar':True,'showXYPlane':False,'showYZPlane':False,'showZXPlane':False})
		self.view.show()


