import numpy as np
import pylab as plt
import bl

 #nxt,nxte,nxs,npt,iswpt
  #169  135   20  101   25 
 #rl,xtr,etae,vgp,deta(1),p2(1)
 #4000000.000   0.00010   8.00000   1.12000   0.01500   1.00000

def amean(ns,nd,x,yy,it):
	y=np.copy(yy)
	d=np.zeros_like(x)
	if ((nd-ns) > 2 and it>0):
		for k in range(0,it):
			for i in range(ns,nd):
				d[i] = y[i]
			for i in range(ns+1,nd-1):
				f1    = (x[i+1] -x[i])/(x[i+1] - x[i-1])
				y[i]  = 0.5*(f1*d[i-1]+d[i]+(1.0-f1)*d[i+1])
	return y

def isfinite(x,k=20):
	if np.isfinite(x) and abs(x)<k:# and x>0.0:
		return True
	else:
		return False

class IblmStrip(object):
	def __init__(self):
		self.sdir='up'
		self.mxp=401
		self.LRef=1.0
		self.nu=1.45e-5
		self.VRef=60.
		self.pts=np.zeros((self.mxp,2),dtype=np.float32)
		self.uebl=np.zeros(self.mxp,dtype=np.float32)
		self.dlspbl=np.zeros(self.mxp,dtype=np.float32)
		self.vnpbl=np.zeros(self.mxp,dtype=np.float32)
		self.uedispbl=np.zeros(self.mxp,dtype=np.float32)
		self.cfsbl=np.zeros(self.mxp,dtype=np.float32)
		self.nptt=200
		self.nlpt=100
		self.iswpt=35
		self.irestart=0
		self.istop=0
		self.sclReynolds=1.
		
		
	def setPts(self,pts):
		self.nptt=pts.shape[0]
		self.pts.fill(0.0)
		self.pts[0:self.nptt]=pts
		
	def rl(self,scl=None):
		if scl is None:
			scl=self.sclReynolds
		return self.VRef*self.LRef/self.nu*scl
	
	def setUebl(self,uebl):
		self.uebl.fill(0.0)
		self.uebl[0:self.nptt]=uebl
	
	def setArray(self,name,array):
		a=getattr(self,name)
		a.fill(0)
		a[0:self.nptt]=array
		
	def setNlpt(self,nlpt):
		self.nlpt=nlpt
		
	def compute(self):
		
		self.ds=np.zeros(len(self.pts))

		for i in range(1,len(self.ds)):
			self.ds[i]=self.ds[i-1]+np.linalg.norm( self.pts[i]-self.pts[i-1])

		self.uebl=np.copy(amean(0,self.nptt-1,self.ds,self.uebl,30))

		
		#blmain(istop,irestart,nxtein,rlin,cdf,xcbl,ycbl,uebl,dlspbl,vnpbl,uedispbl,cfsbl,fmach,iswptin,iang,[nxtin])
		self.istop=0
		nxtin=self.nptt
		nxtein=self.nlpt+1
		rlin=self.rl()
		cdf=0.0#01
		xcbl=self.pts[:,0].astype(np.float32)
		ycbl=self.pts[:,1].astype(np.float32)
		fmach=np.float32(0.)
		iswptin=np.int32(32)
		iang=np.int32(1)
		
		bl.blmain(self.istop,self.irestart,nxtin,nxtein,rlin,cdf,xcbl,ycbl,self.uebl,self.dlspbl,self.vnpbl,self.uedispbl,self.cfsbl,fmach,iswptin,iang)
		self.uevs=np.copy(bl.blrc01.uevs)
		self.bl=bl
		self.istop=np.copy(bl.blinp3.ierror)
		#bl.blc001.nxt=  np.int32(self.nptt)
		#bl.blc001.nxte= np.int32(self.nlpt)
		#bl.blc001.nxs= np.int32(self.nxs)
		#bl.blc001.npt= np.int32(self.npt)
		#bl.blc001.iswpt= np.int32(self.iswpt)

		#bl.blinp3.rl=  np.float32(self.rl()) 
		#bl.blinp3.xtr=np.float32(self.xtr)
		#bl.blinp3.etae=np.float32(self.etae)
		#bl.blinp3.vgp=np.float32(self.vgp)
		#bl.blgrd1.deta[0]=np.float32(self.deta0)
		#bl.blgtyy.p2[0]= np.float32(self.p20)

		#bl.blrc01.xbl=self.pts[:,0].astype(np.float32)
		#bl.blrc01.ybl=self.pts[:,1].astype(np.float32)
		#bl.blrc01.ubl=self.uebl.astype(np.float32)

		#bl.compute()

	def getBlTheta(self):
		return np.copy(bl.blrc01.theta[0:self.nptt])
	def getBlDisplThick(self):
		return np.copy(bl.blrc01.dls[0:self.nptt])
	def getBlowing(self):
		return np.copy(bl.blres1.vnp[0:self.nptt])

	def setBLParam(self,dom1):
		self.sclReynolds=dom1.sclReynolds
		self.compute()
		
		self.uevs=amean(0,self.nptt,self.ds,self.uevs,30)
		self.vnpbl=amean(0,self.nptt,self.ds,self.vnpbl,30)
		self.dlspbl=amean(0,self.nptt,self.ds,self.dlspbl,30)
		
		#blow=self.getBlowing()
		#dthick=self.getBlDisplThick()
		
		relax=0.15
		if self.istop != 0:
			print('-------------------------BL ERROR at '+str( hex(id(self))))
		if self.istop == 0:
			for i in range(0,len(self.panels)):
				
				p=dom1.getPanelByGid(self.panels[i].gid)
				vloc=np.linalg.norm(p.vloc[[0,1]])
				if not p.isWake:
					#p.w=(-self.vnpbl[i+1]*self.VRef)
					#p.displThick=self.dlspbl[i+1]*self.LRef
					#p.uebl=self.uevs[i+1]*self.VRef
					#p.uebl=self.uevs[i+1]*self.VRef
					#p.ue=   self.uevs[i+1]*self.VRef
					##print p.ue
					##p.uebl=p.ue
					#p.dlspbl=self.dlspbl[i+1]*self.LRef
					#p.vnpbl=self.vnpbl[i+1]*self.VRef
					#p.uedispbl=self.uedispbl[i+1]#*self.VRef

					if isfinite(self.uevs[i+1],vloc*2):
						p.blIsSolved=True
						p.uedispbl= relax*(self.uedispbl[i+1]*self.LRef) + (1.-relax)*p.uedispbl
						p.uebl= relax*(self.uevs[i+1]*self.VRef) + (1.-relax)*p.uebl*self.VRef
						p.ue=np.copy(p.uebl)
						p.ue=np.copy(self.uevs[i+1]*self.VRef)
						
						if isfinite(self.dlspbl[i+1],1.):
							p.displThick= relax*(self.dlspbl[i+1]*self.LRef) + (1.-relax)*p.displThick
						if isfinite(self.vnpbl[i+1],10.):
							p.vnpbl= relax*(self.vnpbl[i+1]*self.VRef) + (1.-relax)*p.vnpbl
							p.w=-p.vnpbl
							
						#p.ue=0.
						#print hex(id(p)),hex(id(p.ue)), p.ue/self.VRef, self.uevs[i+1]
						p.cfsbl=self.cfsbl[i+1]
						p.cf=self.cfsbl[i+1]
					
				if p.isWake:
					
					if self.sdir =='up':
						if isfinite(self.uevs[i+1]):
							p.blIsSolved=True
							p.blUeUp=self.uevs[i+1]*self.VRef
							if isfinite(self.dlspbl[i+1],1.):
								p.blDtUp=self.dlspbl[i+1]*self.LRef
							if isfinite(self.vnpbl[i+1],10.):
								p.blVnUp=-self.vnpbl[i+1]*self.VRef
					
					else:
						if isfinite(self.uevs[i+1]):
							p.blIsSolved=True
							p.blUeDown=self.uevs[i+1]*self.VRef
							if isfinite(self.dlspbl[i+1],1.):
								p.blDtDown=self.dlspbl[i+1]*self.LRef
							if isfinite(self.vnpbl[i+1],10.):
								p.blVnDown=-self.vnpbl[i+1]*self.VRef
		
	#def setBLParam(self):
		#self.compute()
		
		#self.uevs=amean(0,self.nptt,self.ds,self.uevs,30)
		#self.vnpbl=amean(0,self.nptt,self.ds,self.vnpbl,30)
		#self.dlspbl=amean(0,self.nptt,self.ds,self.dlspbl,30)
		
		##blow=self.getBlowing()
		##dthick=self.getBlDisplThick()
		#for i in range(0,len(self.panels)):
			#p=self.panels[i]
			#p.w=(-self.vnpbl[i]*self.VRef)
			#p.displThick=self.dlspbl[i]
			#p.uebl=self.uevs[i]*self.VRef
			##p.uebl=self.uebl[i]*self.VRef
			#p.ue=self.uevs[i]*self.VRef
			#p.dlspbl=self.dlspbl[i]
			#p.vnpbl=self.vnpbl[i]*self.VRef
			#p.uedispbl=self.uedispbl[i]#*self.VRef
			#p.cfsbl=self.cfsbl[i]
			#p.cf=self.cfsbl[i]
			#if p.isWake:
				#if self.sdir =='up':
					#p.blDtUp=self.dlspbl[i]
					#p.blUeUp=self.uevs[i]*self.VRef
					#p.blVnUp=-self.vnpbl[i]#*self.VRef
				#else:
					#p.blDtDown=self.dlspbl[i]
					#p.blUeDown=self.uevs[i]*self.VRef
					#p.blVnDown=-self.vnpbl[i]#*self.VRef
	