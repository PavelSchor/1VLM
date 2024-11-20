import sys
sys.path.insert(0, '/usr/local/lib/python2.7/site-packages/')

import numpy as np
np.setbufsize(8192*2)
#from panelMethod.vortexLines import vortexLine ,vortexPolyLine
from vlm.env import GlobalEnvironment
from utils.b3Vect import B3vect 
vect=B3vect()
import  utils.csys as csys

class VlmPanel(GlobalEnvironment):
	def __init__(self,P1,P2,P3,P4):
		self.gid=0
		self.pgid=0
		self.connectivity={'d1':None,'d2':None,'d3':None,'d4':None}
		self.forceVector=np.zeros(3)
		self.dw=np.zeros(3)
		self.G=1.0
		self.P1=P1
		self.P2=P2
		self.P3=P3
		self.P4=P4
		self.INIT_P1=P1.copy()
		self.INIT_P2=P2.copy()
		self.INIT_P3=P3.copy()
		self.INIT_P4=P4.copy()
		self.pts=np.vstack((P1.copy(),P2.copy(),P3.copy(),P4.copy()))
		self.ptsInit=np.vstack((P1.copy(),P2.copy(),P3.copy(),P4.copy()))
		self.csys=np.eye(3)
		self.isWake=False
		self.isLast=False
		self.initGeometry()	

	def initGeometry(self):	
		v12=vect.makeVect(self.pts[0],self.pts[1])
		v13=vect.makeVect(self.pts[0],self.pts[2])
		v14=vect.makeVect(self.pts[0],self.pts[3])
		v23=vect.makeVect(self.pts[1],self.pts[2])
		v24=vect.makeVect(self.pts[1],self.pts[3])
		
		self.boundVector=v12

		p025a=self.pts[0]+vect.setVectSize(v14,0.25*vect.vectorMagn(v14))
		p025b=self.pts[1]+vect.setVectSize(v23,0.25*vect.vectorMagn(v23))

		p075a=self.pts[0]+vect.setVectSize(v14,0.75*vect.vectorMagn(v14))
		p075b=self.pts[1]+vect.setVectSize(v23,0.75*vect.vectorMagn(v23))

		cp=p075a+vect.makeVect(p075a,p075b)*0.5

		gVect=vect.makeVect(p025a, p025b)		
		self.gammaVector=gVect
		#self.dY=np.linalg.norm(self.pts[1]-self.pts[0])
		self.dY=np.linalg.norm(self.gammaVector)
		
		self.controlPoint=self.pts.sum(axis=0)/4.
		self.vortexPoint0=p025a+vect.makeVect(p025a,p025b)*0.5
		self.vortexPoint1=p025a
		self.vortexPoint2=p025b

		self.normalVector=np.cross(v24,v13)/np.linalg.norm(np.cross(v24,v13))
		#/np.linalg.norm(vect.cross(v24,v13))#vect.normal1Vect(v13,v24)
		self.tangentVector=vect.makeVect(p025a,p075a)
		
		self.tangentVector_unit= self.tangentVector/ np.linalg.norm(self.tangentVector)

		self.area=np.linalg.norm( np.cross(self.pts[2]-self.pts[0], self.pts[3]-self.pts[1]) )*0.5
		chord1=self.pts[2]-self.pts[1]
		chord2=self.pts[3]-self.pts[0]
		self.chord=np.mean([chord1,chord2])
		self.gVel=np.zeros(3)
		self.pressureCoef=1.0
#		self.setWake(self.getFreeVelocity())

		self.csys=csys.makeCSYS_from2Vectors(self.tangentVector,self.normalVector,ax=2)

	def getS(self):
		return 0.0
#	def setWake(self,vel):
#		kWake=10e10
#		vel=vel/np.linalg.norm(vel)
#		self.vl0=vortexLine(self.getVortexPoint1(),self.getVortexPoint2())
#
#		self.vl1=vortexPolyLine()
#		self.vl1.addSegment(self.pts[3],self.getVortexPoint1())
#		self.vl1.addSegment(self.pts[3]+vel*kWake,self.pts[3])
#
#		self.vl2=vortexPolyLine()
#		self.vl2.addSegment(self.getVortexPoint2(),self.pts[2])
#		self.vl2.addSegment(self.pts[2],self.pts[2]+vel*kWake)



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
		return self.controlPoint
	
	def getNV(self):
		return self.normalVector

	def setG(self,G):
		self.G=float(G)
		self.vl0.setG(float(G))
		self.vl1.setG(float(G))
		self.vl2.setG(float(G))

	def getG(self):
		return self.G

	def getUVW(self,PT,sf0=1.0,sf1=1.0,sf2=1.0):
		uvw=np.zeros(3)
		uvw=uvw+(self.vl0.getUVW(PT))*sf0
		uvw=uvw+(self.vl1.getUVW(PT))*sf1
		uvw=uvw+(self.vl2.getUVW(PT))*sf2
		
		return uvw

	def getLift(self):
		v= self.getFreeStream()#FreeVelocity()# Stream()
#		vm=vect.vectorMagn(v)#np.linalg.norm(v)
#		vu=v/vm

#		proj=vect.cross(vu,vect.cross(self.boundVector,vu))
#		projm=np.linalg.norm(proj)

		rho=GlobalEnvironment.getRho()
#		boundu=self.boundVector/np.linalg.norm(self.boundVector)
#		self.lu=vect.cross(boundu,vu)
#		self.lu=self.lu/np.linalg.norm(self.lu)
#		projv=np.cross(self.getNV() ,vect.cross(self.getFreeVelocity() , self.getNV()))	
#		if np.dot( self.boundVector, vect.cross(projv,self.getFreeVelocity()- projv) ) < 0:
#		projv=np.cross(self.getNV() ,vect.cross(self.getFreeStream() , self.getNV()))
#		if np.dot( self.boundVector, vect.cross(projv,self.getFreeStream()- projv) ) < 0:
#			sigL=-1.
#		else:
#			sigL=1.
#		sigL=1.
#		print projm
		lift= -rho*np.cross(v,self.boundVector*self.getG())# sigL*self.lu*rho*vm*projm*self.getG()
#		dwm=vect.vectorMagn(self.dw)#np.linalg.norm(self.dw)
#		dwu=self.dw/dwm
#		drag=(vu*rho*dwm*self.getG()*projm)
		return lift#+drag			
		
	def setForce(self):
		rho=GlobalEnvironment.getRho()
#		vel=self.getFreeStream()
		vel=self.getFreeVelocity()
		A=self.getArea()
		self.lift=self.getLift()
		FV=self.getLift()
		force=vect.vectorMagn(FV)#np.linalg.norm(FV)
		cp=(force/A/rho/0.5/np.linalg.norm(self.getFreeVelocity() )**2.0)*np.sign(self.getG())*(-1.)
#		cp=1-(np.linalg.norm(self.getFreeStream())**2. / np.linalg.norm(self.getFreeVelocity())**2.)
		self.setPressureCoef(cp)
		self.forceVector=self.getLift()	
#		self.forceVector=self.getNV()/np.linalg.norm(self.getNV())*cp*self.getArea()
	def getForce(self):
		return self.forceVector

