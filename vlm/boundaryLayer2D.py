import numpy as np
from scipy.interpolate import spline
import scipy.fftpack
from scipy.interpolate import UnivariateSpline
import pylab as plt
def trapz(f,x):
	return np.sum( 0.5*(f[1::]+f[0:-1])*(x[1::]-x[0:-1]))
def thwaites(s,U,Lref,Uref,nu):
	Ue=U/Uref
	n=len(s)
	tht=np.zeros(n)
	theta=np.zeros(n)
	ReTheta=np.zeros(n)
	ReS=np.zeros(n)
	cLambda=np.zeros(n)
	cf=np.zeros(n)
	delta=np.zeros(n)
	H=np.zeros(n)
	CL=np.zeros(n)
	ursum = 0.0
	RE    = Uref*Lref/nu
	f2    = 0.0
	dl=s[1]-s[0]
	dueds=np.gradient(Ue,dl)
	
	
	for i in range(0,n):
		f2=Ue[i]**5
		if i ==0:
			tht[i]= 0.075/dueds[i]
		else:
			ursum=ursum+0.5*(f1+f2)*dl
			const =0.45/(f2*Ue[i])
			tht[i]=const*ursum
		f1 = f2
		theta[i]=np.sqrt(tht[i]/RE)*Lref
		ReTheta[i]=theta[i]*Ue[i]/nu*Uref
		ReS[i]=Ue[i]*s[i]/nu*Uref*Lref
		cLambda[i] = tht[i]*dueds[i]
		if cLambda[i] < 0.0:
			if cLambda[i]+0.107 ==0.0:
				cLambda[i]=-0.1069
			H[i]     = 0.0731/(0.14+cLambda[i])+2.088
			CL[i]    = 0.22+1.402*cLambda[i]+0.018*cLambda[i]/(cLambda[i]+0.107)
		else:
			H[i]     = 2.61-3.75*cLambda[i]+5.24*cLambda[i]**2
			CL[i]    = 0.22+1.57*cLambda[i]-1.8*cLambda[i]**2
		delta[i]=theta[i]*H[i]*Lref
		if i >0:
			cf[i]=2.0*CL[i]/(Ue[i]*theta[i]/Lref*RE)
		if (i > 1) and (cf[i] < 0.0):
			break
		if i >0 and delta[i]<0.0:
			delta[i]=delta[i-1]
			break
		iTrans=i
	return iTrans,theta,delta,H,cf,ReTheta,ReS,tht,cLambda


def fHL(lambdaP):
	H=np.zeros(len(lambdaP))
	for i in range(len(lambdaP)):
		if lambdaP[i]<0.:
			if lambdaP[i] == -0.14:
				lambdaP[i] = -0.139	
			H[i] = 2.088+(0.0731/(lambdaP[i]+0.14))
			
		elif lambdaP[i]>=0:
			H[i] = 2.61-3.75*lambdaP[i]+5.24*lambdaP[i]**2
	return H

def fLL(lambdaP):
	lPlate=np.zeros(len(lambdaP))
	for i in range(len(lambdaP)):
		if lambdaP[i]<0:
			if lambdaP[i]==-0.107:
				lambdaP[i]=-0.106
			lPlate[i] = 0.22+1.402*lambdaP[i] + (0.018*lambdaP[i])/      (lambdaP[i]+0.107)
		if lambdaP[i]>=0:
			lPlate[i] =0.22 + 1.57*lambdaP[i] - 1.8*lambdaP[i]**2
	return lPlate


# entrainment shape factor
def fH1(H):
	#if H <=1.1:
		#return 16.0
	if H<=1.6: 
		return 3.3 + 0.8234*(H-1.1)**(-1.287)
	else:
		return 3.3 + 1.5501*(H-0.6678)**(-3.064)

# shape factor, the inverse function
def fH(H1):
	if H1<=3.32:
		return 3.0
	elif H1>3.32 and H1<5.3:
		return 1.1 + ((H1-3.3)/0.8234)**(-(1./1.287))
	else:
		return 0.6678 + ((H1-3.3)/1.5501)**(-(1./3.064))

# von Karman momentum integral equation    
def F(H1):
	return 0.0306*(H1-3.)**(-0.6169)

def RHS1(theta,H,Ve,nu,dUeds):
	return 0.5*cf(theta,H,Ve,nu) - theta/Ve*(2+H)*dUeds

def RHS2(H1,theta,H,Ve,nu,dUeds):
	return -(H1/theta)*RHS1(theta,H,Ve,nu,dUeds) - (H1/Ve)*dUeds + F(H1)/theta


# skin friction
def cf(theta,H,Ve,nu):
	return 0.246*(10.**((-0.678*H)) )*(theta*Ve/nu)**(-0.268)
	#0.246/((10**(0.678*ht))*((Reo*yo(1)*vo)**0.268))
	#cf = 0.246*(10.^(-0.678*H))*rtheta.^(-0.268);

def differentialOperator(x):
	n=x.shape[0]
	D=np.zeros((n,n))
	d11=x[1]-x[0]
	d12=x[2]-x[0]
	dn1=x[n-2]-x[n-1]
	dn2=x[n-3]-x[n-1]

	D[0,0]= -(d11+d12)/(d11*d12)
	D[0,1]= -(d12)/(d11*(d11-d12))
	D[0,2]=  (d11)/(d12*(d11-d12))

	D[n-1,n-3]= dn1/(dn2*(dn1-dn2))
	D[n-1,n-2]=-dn2/(dn1*(dn1-dn2))
	D[n-1,n-1]=-(dn1+dn2)/(dn1*dn2)
	
	for i in range(1,n-1):
		di1=x[i-1]-x[i]
		di2=x[i+1]-x[i]

		D[i,i-1]= -di2/(di1*(di1-di2))
		D[i,i]  = -(di1+di2)/(di1*di2)
		D[i,i+1]=  di1/(di2*(di1-di2))

	return D

def getDerivative(f,x):
	return np.dot(differentialOperator(x),f)

class BLStrip2D(object):

	def __init__(self):
		self.isReversed=False
		self.panels=[]
		self.wakePoints=np.zeros(2)
		self.mu=1.460e-5
		self.rho=1.225
		self.NI=300
		pass

	def setBlowing(self):
		for i in range(0,len(self.panels)):
			p=self.panels[i]
			p.setBlowingVelocity(self.blowing[i])
			
	def setPoints(self,pts):
		self.pts=pts

	def setDs(self,ds):
		self.ds0=ds-ds[0]
		self.Lref=ds.max()
		self.ds=np.linspace(0,self.ds0[-1],self.NI)#/self.Lref
		
	def setVelocities(self,v):
		#self.Ue0=v
		self.Ue0=np.zeros(len(self.panels))
		for i in range(0,len(self.panels)):
			p=self.panels[i]
			self.Ue0[i]=p.Ue
			
		self.Ue=np.interp(self.ds,self.ds0,self.Ue0)#spline(self.ds0,v,self.ds)
		self.Ue[0]=np.abs(self.Ue[0])
		
	def setVelocitiesFromCP(self,v0,cp):
		self.Ue0=v0
		self.Ue=v0*(1.-cp)**0.5

	def addWakePoints(self,wpts):
		self.wakePoints=wpts

	def exportUe(self,fname):
		np.savetxt(fname,self.Ue0)
	
	def exportDs(self,fname):
		np.savetxt(fname,self.ds0)
		
	def exportPts(self,fname):
		np.savetxt(fname,self.pts[:,[0,2]])
	
	def setDisplThickAnalytic(self):
		self.Vref=np.linalg.norm(self.panels[0].getFreeVelocity())
		deltaOld=np.zeros(len(self.panels))
		#for i in range(0,len(self.panels)):
			#deltaOld[i]=self.panels[i].displThick
		self.deltaOld=np.interp(self.ds,self.ds0,deltaOld)
		#self.Ue/=self.Vref
		#self.ds/=self.Lref
		REE=0.1
		
		rho=self.rho
		mu=self.mu
		self.nu = self.mu/self.rho                 # kinematic viscosity
		nu=self.nu
		self.Re=self.Lref*self.Vref/nu
		
		n=self.ds.shape[0]
		self.H1=np.zeros(n)
		self.dl=np.diff(self.ds)[0]
		self.RE=self.Vref*self.Lref*self.rho/self.mu
		#theta = np.zeros(len(self.ds),dtype=float)

		self.dUeds = np.gradient(self.Ue,self.dl)
		#theta[0] = np.sqrt(0.075*nu/dUeds[0])

		#for i in range(1,len(theta)):
			#theta[i]=np.sqrt( (theta[i-1]**2)*(self.Ue[0]/self.Ue[i])**6 + 0.45*nu/(self.Ue[i]**6)*trapz(self.Ue[0:i+1]**5,self.ds[0:i+1]))
		#theta[2] = theta[1]
		#lambdaP= np.power(theta,2)/nu*np.gradient(self.Ue,self.dl)
		## calculating the parameter l
		#lPlate = fLL(lambdaP)
		## calculating the shape factor from lambda
		#self.H = fHL(lambdaP)

		#self.theta=theta
		#self.dUeds=dUeds
		#self.lambdaP=lambdaP
		#self.lPlate=lPlate

		#cfPlate = 2*lPlate/(self.Re*theta)

		#ReTheta = self.Ue*theta/nu     # Re based on momentum thickness
		#ReS = (rho*self.Ue*self.ds)/nu             # Re based on position

		self.iLTrans,self.theta,self.delta,self.H,self.cf,self.ReTheta,self.ReS,tht,cLambda=thwaites(self.ds,self.Ue,self.Lref,self.Vref,self.nu)
		# ### Michael's Transition Criterion
		# the criterion is entirely based on the reynolds numbers computed below
		mc = 1.174*(np.power(self.ReS,0.46) + 22400*np.power(self.ReS,(-0.54)) )     # transition criterion
		self.mc=mc
		self.iTrans=0
		#self.HGrad=np.gradient(self.H,self.dl)
		for i in range(1,self.iLTrans):
			self.iTrans=i
			self.H1[i]=fH1(self.H[i])
			if self.mc[i]<self.ReTheta[i]:
				print 'Transition point at: ', self.ds[i]
				break
			#if self.HGrad[i] <0. and i>self.iLTrans/2:
				#break
		
		#try:
		#HSep=3.0
		for i in range(self.iTrans-2,len(self.ds)-1):


			self.theta[i+1] = self.theta[i] + self.dl*RHS1(self.theta[i],self.H[i], self.Ue[i],self.nu,self.dUeds[i])
			self.H1[i+1] = self.H1[i] + self.dl*RHS2(self.H1[i],self.theta[i],self.H[i],  self.Ue[i],self.nu,self.dUeds[i])
			self.H[i+1] = fH(self.H1[i+1])
			self.cf[i] = cf(self.theta[i],self.H[i],self.Ue[i],self.nu)
			#if self.H[i]>=HSep and i>len(self.ds)/2:
			if self.H[i+1]<self.H[i] and i>len(self.ds)*0.6:
				dtheta=self.theta[i]-self.theta[i-1]
				self.H[i::]=self.H[i-1]
				for j in range(i,n-1):
					self.theta[j+1]=self.theta[j]+dtheta#.dl*dtheta
					self.cf[j] = cf(self.theta[j],self.H[j],self.Ue[j],self.nu)
				break

		
		self.delta = self.H*self.theta
		try:
			self.delta[np.negative(np.isfinite(self.delta))]=self.delta[np.isfinite(self.delta)].min()
			#self.theta[np.negative(np.isfinite(self.theta))]=0.0#self.theta[np.isfinite(self.theta)].min()
		except:
			print self.delta
		
		damping=1.
		#deltaNew=np.interp(self.ds0,self.ds,self.delta)
		self.deltaDiff=self.delta-self.deltaOld
		self.deltaNew=self.deltaOld+self.deltaDiff*damping
		self.blowing0=np.gradient(self.Ue*self.deltaNew,self.dl)
		self.blowing=np.interp(self.ds0,self.ds,self.blowing0)

		#self.blowing=np.zeros_like(self.ds0)

		delta0=np.interp(self.ds0,self.ds,self.delta)
		for i in range(0,len(self.panels)):
			p=self.panels[i]
			p.displThick=delta0[i]
			
		#deltaNew=np.interp(self.ds0,self.ds,self.delta)
		#deltaDiff=deltaNew-deltaOld
		#self.deltaNew=deltaNew
		
		#if not self.isReversed:
			#for i in range(0,len(self.panels)-1):
				#p=self.panels[i]
				#p.displThick+=deltaDiff[i]
				#p2=self.panels[i+1]
				#p.pts[2]-=p.normalVect*deltaDiff[i]
				#p.pts[3]-=p.normalVect*deltaDiff[i]
				#p2.pts[0]-=p.normalVect*deltaDiff[i]
				#p2.pts[1]-=p.normalVect*deltaDiff[i]

		#if self.isReversed:
			#for i in range(0,len(self.panels)-1):
				#p=self.panels[i]
				#p2=self.panels[i+1]
				#p.displThick+=deltaDiff[i]
				#p.pts[0]-=p.normalVect*deltaDiff[i]
				#p.pts[1]-=p.normalVect*deltaDiff[i]
				#p2.pts[2]-=p.normalVect*deltaDiff[i]
				#p2.pts[3]-=p.normalVect*deltaDiff[i]

		#plt.plot(self.ds,self.deltaNew)