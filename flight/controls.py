import numpy as np
from utils.bCad import BCad

def actualDeflection(inp,lim):
	if np.signbit(inp):
		return abs(inp)*lim[0]
	else:
		return abs(inp)*lim[1]
	
	
class ControlSurfaces(object):

		
	def __init__(self):
		self.bc=BCad()
		self.inpAileron=0.0
		self.inpElevator=0.0
		self.inpRudder=0.0
		self.inpThrottle=0.0
		
		self.dAileronL=0.0
		self.dAileronR=0.0
		self.dElevator=0.0
		self.dRudder=0.0
		
		self.old_dAileronL=0.0
		self.old_dAileronR=0.0
		self.old_dElevator=0.0
		self.old_dRudder=0.0
		
		
		self.limitAileron=(-30.,20.)
		self.limitElevator=(-30.,30.)
		self.limitRudder=(-30.,30.)
	
		self.tAileron=None
		self.tElevator=None
		self.tRudder=None
		self.thrustVector=np.array([-1.,0,0])

		self.limiterAileron=None
		self.limiterElevator=None
		self.limiterRudder=None
		
		self.saveHistory=False
		
		self.hTime=[]
		self.hAileron=[]
		self.hElevator=[]
		self.hRudder=[]
		
	def runScript(self,fname):
		execfile(fname)
	
	def setDeflectionsInTime(self,t):
		l=['Aileron','Elevator','Rudder']
		for i in l:
			if not getattr(self,'t'+i) is None:
				arr=getattr(self,'t'+i)
				deflection=np.interp(t,arr[:,0],arr[:,1])
				try:
					limiter=getattr(self,'limiter'+i)
					deflection=limiter.compute(deflection)
				except:
					pass
				setattr(self,'inp'+i,deflection)
				
		self.deformMesh()

		if self.saveHistory:
			for i in l:
				a=getattr(self,'h'+i)
				a.append(getattr(self,'inp'+i))
				
	def computeDeflections(self):
		self.dAileronL=actualDeflection(self.inpAileron,self.limitAileron)
		self.dAileronR=actualDeflection(-self.inpAileron,self.limitAileron)
		self.dElevator=actualDeflection(self.inpElevator,self.limitElevator)
		self.dRudder=actualDeflection(self.inpRudder,self.limitRudder)
		
	def deformMesh(self):
		self.computeDeflections()
		self.deflectAileronL(self,self.dAileronL-self.old_dAileronL)
		self.deflectAileronR(self,self.dAileronR-self.old_dAileronR)
		self.deflectElevator(self,self.dElevator-self.old_dElevator)
		self.deflectRudder(self,self.dRudder-self.old_dRudder)
		self.old_dAileronL=self.dAileronL
		self.old_dAileronR=self.dAileronR
		self.old_dElevator=self.dElevator
		self.old_dRudder=self.dRudder

	def getDeflectionDict(self):
		return {"dAileronR": self.dAileronR,"dAileronL": self.dAileronL ,"dElevator": self.dElevator, "dRudder": self.dRudder}
	
	def loadControlInputs(self,t,tAileron,tElevator,tRudder):
		
		def processControlInput(t,tInp0,limits):
			tInp=np.copy(tInp0)
			tInp[tInp<0]/=abs(limits[0])
			tInp[tInp>0]/=abs(limits[1])
			return np.hstack(( t[np.newaxis].T, tInp[np.newaxis].T))
		
		#self.tAileron= processControlInput(t,tAileron,self.limitAileron)
		#self.tElevator= processControlInput(t,tElevator,self.limitElevator)
		#self.tRudder= processControlInput(t,tRudder,self.limitRudder)
		if tAileron is not None:
			self.tAileron= processControlInput(t,tAileron,self.limitAileron)
		else:
			print('Not loaded control input : Aileron')
		if tElevator is not None:
			self.tElevator= processControlInput(t,tElevator,self.limitElevator)
		else:
			print('Not loaded control input : Elevator')
		if tRudder is not None:
			self.tRudder= processControlInput(t,tRudder,self.limitRudder)
		else:
			print('Not loaded control input : Rudder')