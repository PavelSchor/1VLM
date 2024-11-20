import numpy as np
from scipy.optimize import curve_fit

import copy
import ctypes
from numpy.ctypeslib import ndpointer

from panelMethod.env import GlobalEnvironment

from utils.b3Vect import B3vect 
vect=B3vect()

#from panelMethod.vortexLines import vortexLine ,vortexPolyLine
from panelMethod.pyDoubletSource import *
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

###	print('panelLib not found!')
############################# vortexLine getUVW ###########################################################################
##testlib.vortexLine_getUVW.restype = ctypes.c_void_p
##testlib.vortexLine_getUVW.argtypes = [ndpointer(ctypes.c_double),ndpointer(ctypes.c_double),ndpointer(ctypes.c_double),ndpointer(ctypes.c_double), ctypes.c_double ]
##testlib.vortexRing_getUVW.restype = ctypes.c_void_p
###(double* uvw, double* P1, double* P2 , double* P3 , double* P4 ,double* PT , double G)
##testlib.vortexRing_getUVW.argtypes = [ndpointer(ctypes.c_double),ndpointer(ctypes.c_double),ndpointer(ctypes.c_double),ndpointer(ctypes.c_double),ndpointer(ctypes.c_double),ndpointer(ctypes.c_double),  ctypes.c_double ]

######################### doubletSourcePanel getPHI ##################################################################
testlib.doubletSourcePanel_getPHI.restype = ctypes.c_double;
testlib.doubletSourcePanel_getPHI.argtypes = [ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),	ctypes.c_int,
		ndpointer(ctypes.c_double),
		ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
		ndpointer(ctypes.c_double),
		ctypes.c_double,
		ctypes.c_double]

######################### doubletPanel getPHI ##################################################################
testlib.doubletPanel_getPHI.restype = ctypes.c_double;
testlib.doubletPanel_getPHI.argtypes = [ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),	ctypes.c_int,
		ndpointer(ctypes.c_double),
		ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
		ndpointer(ctypes.c_double),
		ctypes.c_double]

######################### sourcePanel getPHI ##################################################################
#double sourcePanel_getPHI(double*pts,int n, double*ORIG, double*csysLoc, double*PT,int showPTS){

testlib.sourcePanel_getPHI.restype = ctypes.c_double;
testlib.sourcePanel_getPHI.argtypes = [ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),	ctypes.c_int,
		ndpointer(ctypes.c_double),
		ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
		ndpointer(ctypes.c_double),
		ctypes.c_int]


############################ doubletSourcePanel getUVW ################################################################
#void doubletSourcePanel_getUVW(double*res,double*pts,int n, double*ORIG, double*csysLoc, double*PT,double G,double S,double RC,int*activeEdges)

testlib.doubletSourcePanel_getUVW.restype = ctypes.c_void_p;
testlib.doubletSourcePanel_getUVW.argtypes = [ndpointer(ctypes.c_double),
		ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
		ctypes.c_int,
		ndpointer(ctypes.c_double),
		ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
		ndpointer(ctypes.c_double),
		ctypes.c_double,
		ctypes.c_double,
		ctypes.c_double,
		ndpointer(ctypes.c_int)]
############################ make CSYS from two vectors 12 #################
testlib.vect_makeCSYS_from2Vectors12.restype = ctypes.c_void_p;
testlib.vect_makeCSYS_from2Vectors12.argtypes = [ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),ndpointer(ctypes.c_double),ndpointer(ctypes.c_double)]
#######################################################################

############################ make CSYS from two vectors 13 #################
testlib.vect_makeCSYS_from2Vectors13.restype = ctypes.c_void_p;
testlib.vect_makeCSYS_from2Vectors13.argtypes = [ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),ndpointer(ctypes.c_double),ndpointer(ctypes.c_double)]
#######################################################################

############################ transform vector #########################
#void vect_transformVector(double*res,double*v,double*targetCSYS,double*sourceCSYS){
testlib.vect_transformVector.restype = ctypes.c_void_p;
testlib.vect_transformVector.argtypes = [ndpointer(ctypes.c_double),ndpointer(ctypes.c_double),ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS')]



def uniqueRows(a):
	a = np.ascontiguousarray(a)
	unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
	return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))



class doubletSourcePanel(GlobalEnvironment):
	def __init__(self,P1,P2,P3,P4):
		self.gid=0
		self.pgid=0
		self.vloc=np.zeros(3)
		self.connectivity={'d1':None,'d2':None,'d3':None,'d4':None}
		self.firstPanel=None
		self.lastPanel=None
		self.isWake=False
		self.isFirst=False
		self.isLast=False
		self.wakes=None
		self.stripPanels=[]
		self.forceVector=np.zeros(3)
		self.dw=np.zeros(3)
		self.G=1.0
		self.S=0.0
		self.w=0.0
		self.cf=0.0
		self.dP=0.0
		self.jumpPot=0.0
		self.displThick=0.0
		self.Ue=0.0
		self.pts=np.vstack((P1,P2,P3,P4))
		self.P1=self.pts[0]
		self.P2=self.pts[1]
		self.P3=self.pts[2]
		self.P4=self.pts[3]	
		self.initGeometry()	
		self.activeEdges=np.ones(4,dtype=np.int32)
		self.ue=0.
		self.uebl=0.
		self.dlspbl=0.
		self.vnpbl=0.0
		self.uedispbl=0.
		self.cfsbl=0.
		self.blDtUp=0.
		self.blDtDown=0.
		self.blUeUp=0.
		self.blUeDown=0.
		self.blVnUp=0.0
		self.blVnDown=0.
		self.blIsSolved=True
		self.blWakeSolveVelocity=False
		self.swapGamaGradients=False
		self.vRel=np.zeros(3)
		self.vFarWake=np.zeros(3)	
	def initGeometry(self):	

	
#		self.pts=uniqueRows(self.pts)
		self.pts.astype(np.float64)
		v12=vect.makeVect(self.P1,self.P2)
		v13=vect.makeVect(self.P1,self.P3)
		v14=vect.makeVect(self.P1,self.P4)
		v23=vect.makeVect(self.P2,self.P3)
		v24=vect.makeVect(self.P2,self.P4)
		v34=vect.makeVect(self.P3,self.P4)
		self.boundVector=v12.astype(np.float64)

		p025a=self.P1+vect.setVectSize(v14,0.5*vect.vectorMagn(v14))
		p025b=self.P2+vect.setVectSize(v23,0.5*vect.vectorMagn(v23))

		p050a=self.P1+vect.setVectSize(v14,0.5*vect.vectorMagn(v14))
		p050b=self.P2+vect.setVectSize(v23,0.5*vect.vectorMagn(v23))

		p075a=self.P1+vect.setVectSize(v14,0.75*vect.vectorMagn(v14))
		p075b=self.P2+vect.setVectSize(v23,0.75*vect.vectorMagn(v23))


		p050le=self.P1+vect.setVectSize(v12,0.5*vect.vectorMagn(v12))
		p050te=self.P3+vect.setVectSize(v34,0.5*vect.vectorMagn(v34))

		self.vortexPoint0=self.P1#p025a+vect.makeVect(p025a,p025b)*0.5
		self.vortexPoint1=self.P2#p025a
		self.vortexPoint2=self.P3#p025b
		self.vortexPoint3=self.P4
		

		self.vx=p050le-p050te
		self.vy=p050a-p050b
		
		self.vxHalf=self.vx*0.5
		self.vyHalf=self.vy*0.5
		self.normalVect=vect.cross(v24,v13)#/np.linalg.norm(vect.cross(v24,v13))#vect.normal1Vect(v13,v24)
		self.normalVect/=np.linalg.norm(self.normalVect)
		self.tangentVect0=vect.makeVect(p050le,p050te)
		self.tangentVect=self.tangentVect0/np.linalg.norm(self.tangentVect0)
		self.tangentVect=self.tangentVect.astype(np.float64)

		uniquePts=uniqueRows(self.pts)
		if np.any(np.isnan(self.tangentVect)):
			self.tangentVect0=vect.makeVect(uniquePts[0],uniquePts[len(uniquePts)-1])
			self.tangentVect=self.tangentVect0/np.linalg.norm(self.tangentVect0)
			self.tangentVect=self.tangentVect.astype(np.float64)
			self.vx=-self.tangentVect0
			self.vxHalf=self.vx*0.5
		self.length=np.linalg.norm(self.tangentVect0)

		cp=np.sum(uniquePts,axis=0)/len(uniquePts)
#		cp=(self.P1+self.P2+self.P3+self.P4)/4.
#		cp+=(self.normalVect/np.linalg.norm(self.normalVect)*1e-8)
		self.controlPoint=cp.astype(np.float64)
#		self.midPoint=p050a+vect.makeVect(p050a,p050b)*0.5
		self.midPoint=cp#(self.P1+self.P2+self.P3+self.P4)/4.
		gVect=vect.makeVect(p025a, p025b)		
		self.gammaVector=gVect


		self.area=0.5*vect.vectorMagn(vect.cross(v13,v24))

		self.csys=np.zeros((3,3),dtype=np.float64)

		testlib.vect_makeCSYS_from2Vectors13(self.csys,self.tangentVect,self.normalVect)
		self.gVel=np.zeros(3)
		self.pressureCoef=0.0
		if np.dot(self.getFreeVelocity(),self.tangentVect)<0.:
			self.dragVect=self.tangentVect*-1.
		else:
			self.dragVect=self.tangentVect

		######### pyDoubletSource
		#self.pds=pyDSPanel()
		#self.pds.setOrig(self.midPoint)
		#self.pds.setCSYS(self.csys)
		#self.pds.setPoints(self.pts)
		
	
	def transform2loc(self,orig,pts):
#		ptsl=np.zeros(pts.shape)
		ptsl=pts-orig
		if len(ptsl.shape) >1:
			for i in range(0,ptsl.shape[0]):
				ptsl[i]=np.dot(self.csys,ptsl[i])
		if len(ptsl.shape) == 1:
			ptsl=np.dot(self.csys,ptsl)
		return ptsl

	def transform2glob(self,orig,pts):
#		ptsl=np.zeros(pts.shape)
		ptsl=pts-orig
		if len(ptsl.shape) >1:
			for i in range(0,ptsl.shape[0]):
				ptsl[i]=np.dot(np.linalg.inv(self.csys),ptsl[i])
		if len(ptsl.shape) == 1:
			ptsl=np.dot(np.linalg.inv(self.csys),ptsl)
		return ptsl


	def getVORING_UVW(self,PT):
		return 0#self.vl1.getUVW(PT)*self.G

	def connectFlowDomain(self,dom):
		self.dom=dom

	def addWakePanel(self,wPanel):
		self.hasWake=wPanel.gid
		
	def setFreeStream(self,vel):
		self.gVel=vel

	def getArea(self):
		return self.area

	def getGammaVector(self):
		return self.gammaVector

	def setPressureCoef(self,cp):
		self.pressureCoef=cp

	def getPressureCoef(self):
		return self.pressureCoef

	def getFreeStream(self):
		return self.gVel
	
	def getVortexPoint0(self):
		return self.vortexPoint0


	def getVortexPoint1(self):
		return self.vortexPoint1

	def getVortexPoint2(self):
		return self.vortexPoint2

	def getCpoint(self):
		#if self.isFirst or self.isLast:
			#return self.controlPoint-self.getNV()*self.displThick
		#else:
		return self.controlPoint
	def getBlCpoint(self,frac=1.):	
		return self.controlPoint-self.getNV()*self.displThick*frac
	
	def getNV(self):
		return self.normalVect

	def setG(self,G):
		self.G=np.copy(float(G))
	
	def getG(self):
		return self.G

	def setS(self,S):
		self.S=np.copy(float(S))
	
	def getS(self):
		return self.S

	def setBlowingVelocity(self,w):
		self.w=w

	def getBlowingVelocity(self):
		if not self.isWake:
			return self.w
		else:
			return 	self.blVnUp - self.blVnDown

	def getUVW(self,PT,RC=0.0):
		uvw=np.zeros(3,dtype=np.float64)
		testlib.doubletSourcePanel_getUVW(uvw,self.pts,self.pts.shape[0],self.controlPoint,self.csys,PT.astype(np.float64),np.float64(-self.getG()),np.float64(self.getS()),np.float64(RC),self.activeEdges)
		return uvw
		#return self.getDoubletUVW(PT,self.G)+self.getSourceUVW(PT,self.S)

	def getDoubletUVW(self,PT,G,RC=0.0):
		uvw=np.zeros(3,dtype=np.float64)
		testlib.doubletSourcePanel_getUVW(uvw,self.pts,self.pts.shape[0],self.controlPoint,self.csys,PT.astype(np.float64),np.float64(G),np.float64(0.0),np.float64(RC))
		return uvw

	def getSourceUVW(self,PT,S):
		uvw=np.zeros(3,dtype=np.float64)
		testlib.doubletSourcePanel_getUVW(uvw,self.pts,self.pts.shape[0],self.controlPoint,self.csys,PT.astype(np.float64),0.0,np.float64(S))
		return uvw

	def getPHI(self,PT):
		return testlib.doubletSourcePanel_getPHI(self.pts,self.pts.shape[0],self.controlPoint,self.csys,PT.astype(np.float64),np.float64(self.G), np.float64(self.S));

	def getDoubletPHI(self,PT,G):
		return testlib.doubletSourcePanel_getPHI(self.pts,self.pts.shape[0],self.controlPoint,self.csys,PT.astype(np.float64),np.float64(G),np.float64(0.0))

	def getSourcePHI(self,PT,S,showPTS=0):
		return testlib.doubletSourcePanel_getPHI(self.pts,self.pts.shape[0],self.controlPoint,self.csys,PT.astype(np.float64),np.float64(0.0),np.float64(S),showPTS);
		#return testlib.sourcePanel_getPHI(self.pts,self.pts.shape[0],self.controlPoint,self.csys,PT.astype(np.float64),showPTS);
	def getLift(self):
		v= self.getFreeStream()#FreeVelocity()# Stream()
		rho=GlobalEnvironment.getRho()
		G=self.getG()
		if self.connectivity['d3']:
			p=self.dom.getPanelByGid(self.connectivity['d3'])
			G=G-p.getG()
		if self.hasWake2:
			p=self.dom.getPanelByGid(self.hasWake2)
			G=p.getG()-G
		lift= -rho*np.cross(v,self.boundVector*G)# sigL*self.lu*rho*vm*projm*self.getG()
		return lift#+drag			

	def _getDoubletDerivative(self,p1,p2):
		dG=p1.getG()-p2.getG()
		pt1l=self.transform2loc(self.midPoint,p1.midPoint)
		pt2l=self.transform2loc(self.midPoint,p2.midPoint)
		vv=pt2l-pt1l
		vu=vv/np.linalg.norm(vv)
		vr=vu*(dG/np.linalg.norm(vv))
		return vr

	def setForce(self):
		def evalGamaGradient(self,kk):
			dir1='vxHalf'
			dir2='vyHalf'
			if self.swapGamaGradients:
				dir2='vxHalf'
				dir1='vyHalf'	
			isPanel=np.zeros(3,dtype=np.bool)
			dg0=np.zeros(3)
			ds0=np.zeros(3)
			isPanel[1]=True
			dg0[1]=self.getG()
			for k in kk:
				if  self.connectivity[k] is not None:
					p1=self.dom.getPanelByGid(self.connectivity[k])
					isPanel[kk[k]]=True
					dg0[kk[k]]=p1.getG()

					if k == 'd1' or k == 'd3':
						ds0[kk[k]]=np.linalg.norm(getattr(self,dir1))+np.linalg.norm(getattr(p1,dir1))
					if k == 'd2' or k == 'd4':
						ds0[kk[k]]=np.linalg.norm(getattr(self,dir2))+np.linalg.norm(getattr(p1,dir2))
					if (k == 'd2') or  (k == 'd1'):
							ds0[kk[k]]*=-1.
			ds=np.zeros(np.sum(isPanel))
			dg=np.zeros(np.sum(isPanel))
			j=0
			for i in range(0,3):
				if isPanel[i]:
					ds[j]=ds0[i]
					dg[j]=dg0[i]
					j+=1
			A = np.vstack([ds, np.ones(len(ds))]).T
			m, c = np.linalg.lstsq(A, dg)[0]
			
			if np.sum(isPanel)==3:
				ds+=ds[0]
				def fit_func(x, a, b,c):
					return a*x**2 + b*x+c
				
				params= curve_fit(fit_func, ds,dg)
				[a, b, c] = params[0]
				m=2*a*ds[1]+b
				
			
			return m
		
		if not (self.isWake):# or self.isFirst or self.isLast):

			freeVelLoc=self.transform2loc(np.zeros(3),self.getFreeVelocity())
			self.fvl=copy.copy(freeVelLoc)
			freeVelLoc[2]-=self.S###+self.getBlowingVelocity())
			self.vloc.fill(0.0)
			self.vloc+=freeVelLoc
####			try:
####				if self.isLast or self.isFirst:
####					if self.connectivity['d1'] is not None:
####						p3=self.dom.getPanelByGid(self.connectivity['d1'])
####					if self.connectivity['d3'] is not None:
####						p3=self.dom.getPanelByGid(self.connectivity['d3'])
####					dG=self.getG()-p3.getG()
####					dl=np.linalg.norm(p3.vxHalf)+np.linalg.norm(self.vxHalf)
####					
####					self.dG=dG
####					self.dl=dl
####					#### UGLY hack, will probably fail.. the forward/backward difference must be evaluated properly!!
####					if self.isLast:
####						self.vloc[0]+=np.abs(dG/dl)*np.sign(self.vloc[0])
####					if self.isFirst:
####						self.vloc[0]+=np.abs(dG/dl)*np.sign(self.vloc[0])*-1
####			except:
####				self.vloc[0]=np.linalg.norm(self.getFreeVelocity())
####				print('Panel: '+str(self.gid)+' trailing edge connectivity failed')
####				
####				
####			if not (self.isLast or self.isFirst):
####				if self.connectivity['d1'] is not None and self.connectivity['d3'] is not None:
####					
####					p1=self.dom.getPanelByGid(self.connectivity['d1'])
####					p3=self.dom.getPanelByGid(self.connectivity['d3'])
####					dG=p3.getG()-p1.getG()
####					dl=np.linalg.norm(p1.vxHalf)+np.linalg.norm(self.vx)+np.linalg.norm(p3.vxHalf)
####					
####					self.vloc[0]+=dG/dl
####				if self.connectivity['d2'] is not None and self.connectivity['d4'] is not None:
####					try:
####						p2=self.dom.getPanelByGid(self.connectivity['d2'])
####						p4=self.dom.getPanelByGid(self.connectivity['d4'])
####						dG=p4.getG()-p2.getG()
####						dl=np.linalg.norm(p2.vyHalf)+np.linalg.norm(self.vy)+np.linalg.norm(p4.vyHalf)
####						self.vloc[1]+=dG/dl
####					except:
####						pass
			try:
				self.vloc[0]+=evalGamaGradient(self,{'d1':0,'d3':2})
			except:
				pass
			try:
				self.vloc[1]+=evalGamaGradient(self,{'d2':0,'d4':2})
			except:
				pass
			cp=1-((np.linalg.norm(self.vloc))**2. / np.linalg.norm(self.getFreeVelocity())**2.)
			if self.blIter>0 and self.displThick>0.0 and self.blIsSolved:
				cp=1-((self.ue   )**2. / np.linalg.norm(self.getFreeVelocity())**2.)
	
			self.setFreeStream(self.transform2glob(np.zeros(3),self.vloc))
			self.setPressureCoef(cp)
			self.forceVector=self.getNV()*cp*(self.getArea()*0.5*self.getRho()*np.linalg.norm(self.getFreeVelocity())**2.)
			#self.forceVector+=self.dragVect*self.cf*(self.getArea()*0.5*self.getRho()*np.linalg.norm(self.getFreeVelocity())**2.)
			sgn=np.sign(np.dot(self.getForce(),self.getNV()))
			self.dP=sgn*np.linalg.norm(self.getForce())/self.getArea()
		else:
			self.setPressureCoef(0.0)
	
	def getForce(self):
		if self.isWake:
			return np.zeros(3)
		else:
			return self.forceVector
	
	def getForceAtPt(self,PT):
		f=np.zeros(6)
		f[0:3]=self.getForce()
		r=self.midPoint -PT
		f[3:6]=np.cross(r,self.getForce())
		return f

	
	def getPressureLoad(self):
		return self.getPressureCoef()*0.5*self.getRho()*np.linalg.norm(self.getFreeVelocity())**2.

	def getSubCentroids(self,m,n):
		p0=np.zeros((m,3))
		p1=np.zeros((m,3))
		v2=np.zeros((m,3))
		res=np.zeros((m*n,3))
		uniquePts=uniqueRows(self.pts)
		v0=uniquePts[1]-uniquePts[0]
		v1=uniquePts[-1]-uniquePts[2]
		l0=np.linspace(-0,1.0,m+2)
		l1=np.linspace(-0,1.0,n+2)
		for i in range(0,m):
			p0[i]=uniquePts[0]+v0*l0[i+1]
			p1[i]=uniquePts[2]+v1*l0[i+1]
			v2[i]=p1[i]-p0[i]
		for i in range(0,m):
			for j in range(0,n):
				res[i*n+j]=p0[i]+v2[i]*l1[j+1]
		return res