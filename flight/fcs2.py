import numpy as np
from scipy.interpolate import interp1d                                                                                                                                                 

class ControlInputs(object):
	def __init__(self):
		self.aileron=np.array([0.])
		self.elevator=np.array([0.])
		self.rudder=np.array([0.])
		self.throttle=np.array([0.])
		self.flap=np.array([0.])
	def initAll(self):
		self.aileron.fill(0)
		self.elevator.fill(0)
		self.rudder.fill(0)
		self.throttle.fill(0)
		self.flap.fill(0)


class ControlDeflection(object):
	def __init__(self):
		self.rng=np.array((-30.,30.))
		self.inp=0.0
		self.computeRng()
		self.cmds=[]

	def computeRng(self):
		self.mid=self.rng.mean()
		# self.delta=self.rng.max()-self.rng.min()
		self.delta=self.rng[1]-self.rng[0]
		self.dh=self.delta/2.

	def setRng(self,rng):
		self.rng[0]=rng[0];
		self.rng[1]=rng[1];
		self.computeRng()

	def getDeflection(self,inp=None):
		if inp==None:
			inp=getattr(self.inputs,self.cmds[0])
		d=self.mid+self.dh*inp
		if type(d)==np.ndarray:
			d=d[0]
		return d

class CompoundControlDeflection(object):
	def __init__(self):
		self.surfaces=[]
		self.cmds=[]
		self.dmin=0.
		self.dmax=0.

	def getDeflection(self):
		d=0.
		for i in range(0,len(self.surfaces)):
			ci=getattr(self.inputs,self.cmds[i])
			d+=self.surfaces[i].getDeflection(ci)
		if d<self.dmin:
			d=self.dmin
		if d>self.dmax:
			d=self.dmax
		if type(d)==np.ndarray:
			d=d[0]
		return d



#
# c1=ControlDeflection()
# c1.setRng([-30,30])
# c2=ControlDeflection()
# c2.setRng([30,-30])
# c2.getDeflection(-1)
# c2.getDeflection(1)
#
# cc=CompoundControlDeflection()
# cc.dmin=-30.;cc.dmax=30.
# cc.surfaces.append(c1)
# cc.surfaces.append(c2)
#
# ci=ControlInputs()
# cc.inputs=ci
# cc.cmds.append('elevator')
# cc.cmds.append('aileron')
# cc.getDeflection()
#
# ci.aileron=-0.5
# ci.elevator=0.5
# cc.getDeflection()










class FlightControlSystem2(object):
	def __init__(self):
		self.inputs=ControlInputs()
		self.groupDeflections={}
		self.groupSurfaces={}

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
	
	def deflectRudders(self):#,dAileron,dElevator,dRudder):
		for k in self.groupDeflections:
			self.groupSurfaces[k].deflectControlSurface(self.groupDeflections[k].getDeflection())
	
	def setDeflectionsInTime(self,t):
		dAileron,dElevator,dRudder=self.getDeflections1(t)
		self.ruddersMoved=False
		if abs(self.dAileron - dAileron) > 0 or abs(self.dElevator - dElevator) > 0 or abs(self.dRudder - dRudder) > 0:
			self.ruddersMoved =True
		self.dAileron,self.dElevator,self.dRudder=self.getDeflections1(t)
		self.deflectRudders(self.dAileron,self.dElevator,self.dRudder)
	
	def getDeflectionDict(self):
		return {'dAileron':self.dAileron,'dElevator':self.dElevator,'dRudder':self.dRudder}
