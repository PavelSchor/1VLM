import numpy as np

class BeamCrossSection(object):
	def __init__(self):
		self.B=1.
		self.H=1.
		self.b=0.
		self.h=0.
		pass

	def setBHbh(self,B,H,b,h):
		self.B=B; self.H=H; self.b=b; self.h=h;

	def getArea(self):
		return (self.B*self.H) - (self.b*self.h)

	def getIy(self):
		return ((self.B**3  * self.H) - (self.b**3 *self.h))/12.

	def getIz(self):
		return ((self.H**3  * self.B) - (self.h**3 *self.b))/12.
	
	def getSigmaBendingZZ(self,M):
		z=self.H/2.
		return M/self.getIz()*z
		
	def findFlangeHeightZZ(self,M,sigma,hmin=0.01):
		maxIt=100
		tol=1e-3
		z=self.H/2.
		a=0.0
		b=self.H-hmin
		self.h=0.0
		smin=self.getSigmaBendingZZ(M)
		self.h=b
		smax=self.getSigmaBendingZZ(M)
		if sigma >= smax:
			self.h=0.0
			return 
		err=sigma-smax
		sa=smin-sigma
		sb=smax-sigma
		
		i=0
		while (abs(err) > abs(tol)) and i<maxIt:
			i+=1
			c=(a+b)/2.
			self.h=c
			sc=self.getSigmaBendingZZ(M)-sigma
			
			self.h=a
			if np.sign(sc) == np.sign(self.getSigmaBendingZZ(M)-sigma):
				a=c
			else:
				b=c
			err=sc
			self.h=c
		return