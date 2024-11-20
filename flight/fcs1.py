import numpy as np
from scipy.interpolate import interp1d                                                                                                                                                 


class FlightControlSystem1(object):
	def __init__(self):
		self.ruddersMoved=True
		self.isActive_Aileron=True
		self.isActive_Elevator=True
		self.isActive_Rudder=True
		pass
	
	def loadTimeCmds(self,t,tAileron,tElevator,tRudder):
		self.t=t.copy()
		self.tAileron=tAileron
		self.tElevator=tElevator
		self.tRudder=tRudder
		
		self.fAileron=interp1d(self.t,self.tAileron)
		self.fElevator=interp1d(self.t,self.tElevator)
		self.fRudder=interp1d(self.t,self.tRudder)
		
		self.dAileron,self.dElevator,self.dRudder=self.getDeflections1(self.t[0])
		
	def getDeflections1(self,tt):
		dAileron=self.fAileron(tt).flatten()[0]
		dElevator=self.fElevator(tt).flatten()[0]
		dRudder=self.fRudder(tt).flatten()[0]
		return dAileron,dElevator,dRudder
	
	def deflectRudders(self,dAileron,dElevator,dRudder):
		if self.isActive_Aileron:
			self.GRP_Aileron.deflectControlSurface(ctrMask=self.msk_AIL,ax=self.AX_AIL,orig=self.ORG_AIL,d=dAileron)
		if self.isActive_Elevator:
			for i in self.GRP_Elevator:
				i.deflectControlSurface(dElevator)
		if self.isActive_Rudder:
			for i in self.GRP_Rudder:
				i.deflectControlSurface(dRudder)
	
	def setDeflectionsInTime(self,t):
		dAileron,dElevator,dRudder=self.getDeflections1(t)
		self.ruddersMoved=False
		if abs(self.dAileron - dAileron) > 0 or abs(self.dElevator - dElevator) > 0 or abs(self.dRudder - dRudder) > 0:
			self.ruddersMoved =True
		self.dAileron,self.dElevator,self.dRudder=self.getDeflections1(t)
		self.deflectRudders(self.dAileron,self.dElevator,self.dRudder)
	
	def getDeflectionDict(self):
		return {'dAileron':self.dAileron,'dElevator':self.dElevator,'dRudder':self.dRudder}
