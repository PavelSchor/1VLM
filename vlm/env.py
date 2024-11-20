import numpy as np
class GlobalEnvironment(object):
	rho=1.225
	pressure=10132.5
	mu=1.460e-5
	freeVel=np.zeros(3)
	blIter=0

	@classmethod	
	def getRho(cls):
		return cls.rho
	@classmethod
	def getPressure(cls):
		return cls.pressure
	@classmethod
	def getFreeVelocity(cls):
		return GlobalEnvironment.freeVel
	@classmethod	
	def getMu(cls):
		return cls.mu
	@classmethod
	def setRho(cls,rho):
		cls.rho=rho
	@classmethod
	def setPressure(cls,pressure):
		cls.pressure=pressure
	@classmethod
	def setFreeVelocity(cls,vel):
		GlobalEnvironment.freeVel=vel
	@classmethod	
	def setMu(cls,mu):
		GlobalEnvironment.mu=mu	
