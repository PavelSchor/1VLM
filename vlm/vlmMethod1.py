import sys
sys.path.insert(0, '/usr/local/lib/python2.7/site-packages/')


from vlm.vlmPanels import VlmPanel
from numba import jit,prange

import numpy as np


from collections import OrderedDict
import copy
from utils.transform import EulerMatrix
from utils.b3Vect import B3vect 
from utils.airfoils import AeroSurfFactory
vect=B3vect()
from utils.bCad import BCad
bc= BCad()
from vlm.viewer import *

from time import gmtime, strftime
from scipy.sparse.linalg import lgmres

#import affinity
import multiprocessing

from scipy import optimize
#def sig_handler(signum, frame):
	#print('Segmentation Fault')

#signal.signal(signal.SIGSEGV, sig_handler) 

#affinity.set_process_affinity_mask(0,2**multiprocessing.cpu_count()-1)


def projectOnPlane(u,n):
	return u - (np.dot(u,n)/(np.linalg.norm(n)**2.)) *n

def pointInArray(pts,x,tol=1e-5):
	try:
		return np.less(abs(pts-x),np.ones(3)*tol).all(axis=1).any()
	except:
		return False
	
def commonEdgeTwoPanels(A,B):
	nrows, ncols = A.pts.shape
	dtype={'names':['f{}'.format(i) for i in range(ncols)], 'formats':ncols * [A.pts.dtype]}
	C = np.intersect1d(A.pts.view(dtype), B.pts.view(dtype))
	if len(C) < 2:
		return False
	else:
		return True

@jit(nopython=True,nogil=True,parallel=True)
def assembleMatrixA_numba(AIB,PTS,VORTEXPTS,NORMALS):
	m=len(AIB)
	AIB.fill(0.0)
	for i in prange(0,m):
		for j in prange(0,m):
			#aa=0.0
			for k in range(0,4):
				a=PTS[i]-VORTEXPTS[j*4+k]
				b=PTS[i]-VORTEXPTS[j*4+((k+1)%4)]
				an=np.linalg.norm(a)
				bn=np.linalg.norm(b)
				#K= (A+B)/(A*B*(A*B+vect_dot(a,b) ))/4.0/M_PI;
				K=(an+bn)/(an*bn*( an*bn+np.dot(a,b) ))/4./np.pi # np.sum(np.multiply(a,b))))
				V=np.cross(a,b) *  K#[np.newaxis].T
				AIB[i][j]+=  np.sum( np.multiply( V, NORMALS[i]))
				#aa+=  np.sum( np.multiply( V, NORMALS[i]))

@jit(nopython=True,nogil=True,parallel=True)
def computeWakeContribution_numba(AIW,PTS,WAKEPTS,WAKEMASK,NORMALS):
	AIW.fill(0.0)
	m=AIW.shape[0]
	mw=WAKEPTS.shape[1]
	for i in prange(0,m):
		P=PTS[i]
		for j in prange(0,m):
			if WAKEMASK[j]:
				for k in prange(0,mw):
					for l in prange(0,4):
						a=P-WAKEPTS[j][k][l]
						b=P-WAKEPTS[j][k][(l+1)%4]
						an=np.linalg.norm(a)
						bn=np.linalg.norm(b)
						K=(an+bn)/(an*bn*( an*bn+np.dot(a,b) ))/4./np.pi # np.sum(np.multiply(a,b))))
						V=np.cross(a,b) *  K#[np.newaxis].T
						AIW[i][j]+=  np.sum( np.multiply( V, NORMALS[i]))
						
						
class SubSection(object):
	def __init__(self):
		self.nSpanwise=0
		self.nPanChordwise=0
		self.nPanels=0
		self.nWakePanels=0
		self.nWakesChordwise=0
		self.panels=OrderedDict()
		self.wakes=OrderedDict()
		self.refPT=np.zeros(3)

	def addPanel(self,panel):
		self.panels[self.nPanels]=panel
		self.nPanels=self.nPanels+1


	def computePlainGridParams(self):
		pass
		#n=self.nLPanels
		#p0=self.panels[0]
		#for i in range(2,n):
			#if commonEdgeTwoPanels(p0,self.panels[i]):
				#break
		#self.nLPanChordwise=i+1
		#if self.nLPanels%self.nLPanChordwise !=0:
			#p0=self.panels[n-1]
			#j=3
			#for i in range(n-3,-1,-1):
				#if commonEdgeTwoPanels(p0,self.panels[i]):
					#break
				#else:
					#j+=1
			#self.nLPanChordwise=j#+1			
		#self.nSpanwise=self.nLPanels/self.nLPanChordwise
	
		
	def setPGIDS(self):
		for i in range(0,len(self.panels)):
			self.panels[i].pgid=i
	
	def getPanelByPGID(self,i):
		return self.panels[i]


	def getChord(self):
		p0=self.getLiftPanel(0).getCpoint()
		dv=np.zeros(self.getNumberOfLiftPanels())
		for i in range(0,self.getNumberOfLiftPanels()):
			dv[i]=np.linalg.norm(self.getLiftPanel(i).getCpoint()-p0)
		return dv.max()
	


	def getSpanwiseStrips(self,liftOnly=False):

		g=self
		self.nPanChordwise=g.nLPanChordwise+(g.nPanels-g.nLPanChordwise*g.nSpanwise)/g.nSpanwise
		nsec=self.nSpanwise
		npan=self.nPanChordwise

		pg={}
		for s in range(0,nsec):
			pg[s]=PanelGroup()
			pg[s].nPanChordwise=self.nPanChordwise
			pg[s].nLPanChordwise=self.nLPanChordwise
			pg[s].nSpanwise=1
			for i in range(0,npan):
				pg[s].addPanel(self.getPanel(s*npan + i))
				
		if liftOnly:
			npan=self.nLPanChordwise	
			pg={}
			for s in range(0,nsec):
				pg[s]=PanelGroup()
				for i in range(0,npan):
					pg[s].addPanel(self.getLiftPanel(s*npan + i))

		return pg


	
	def setRefPT(self,PT):
		self.refPT=PT

	def getPanelList(self):
		return self.panels.keys()

	def getNumberOfPanels(self):
		return len(self.panels)

	def getPanel(self,i):
		return self.panels[i]

	
	def translate(self,v):
		for i in self.getPanelList():
			self.getPanel(i).pts+=v
			

	def getForceAtPT(self,PT):
		f=np.zeros(3)
		m=np.zeros(3)
		for i in self.getPanelList():
			panelI=self.getPanel(i)
			if not panelI.isWake:
				f=f+panelI.getForce()
				r=panelI.getCpoint() -PT
				m=m+np.cross(r,panelI.getForce())
		return np.hstack((f,m))

	def getForceAtRefPT(self):
		return self.getForceAtPT(self.refPT)

	def getArea(self):
		A=0.0
		for i in self.getPanelList():
			A+=self.getPanel(i).getArea()
		return A		


	def getForceAtPTByGidList(self,PT,l):
		f=np.zeros(3)
		m=np.zeros(3)
		for i in self.getPanelList():
			panelI=self.getPanel(i)
			if panelI.gid in l:
				f=f+panelI.getForce()
				r=panelI.midPoint -PT
				m=m+np.cross(r,panelI.getForce())
		return np.hstack((f,m))

	
	def getPTS2Array(self):
		r=np.zeros((self.nLPanels*4,3))
		for i in range(0,self.nLPanels):
			p=self.getLiftPanel(i)
			for j in range(0,4):
				r[i*4+j]=p.pts[j]
		return r
	
	def deflectControlSurface(self,d):#,ctrMask,ax,orig,d):
		for i in range(0,len(self.ctrMask)):
			if self.ctrMask[i]:
				p=self.panels[i]
				p.pts=bc.rotAx(p.ptsInit,self.AX,self.AX_O,d)
				if p.isLast:
					p.wakes[0].pts[0:2]=p.pts[2:4][::-1].copy()
	def getPGByMask(self,msk):
		pg=PanelGroup()
		for i in range(0,len(self.panels)):
			if msk[i]:
				pg.addPanel(self.panels[i])
		return pg


class PanelGroup(SubSection):
	def addSubsection(self,r,join=True):
		for i in r.getPanelList():
			self.addPanel(r.getPanel(i))
		if join:
			self.nSpanwise+=r.nSpanwise
			self.nLPanChordwise=r.nLPanChordwise
		
	def addPanelGroup(self,pg):
		for i in pg.getPanelList():
			self.addPanel(pg.getPanel(i))


class VLMDomain(object):
	def __init__(self):
		self.fcs=None
		self.VTOT=np.zeros(3)
		self.VEL=np.zeros(3)
		self.COR=np.zeros(3)
		self.OMEGA=np.zeros(3)
		self.nRuns=0
		self.view=PanelViewer3D()
		self.refPT=np.zeros(3)
		self.regions=OrderedDict()
		self.gids=OrderedDict()
		self.nReg=0
		self.nPan=0
		self.nLPan=0
		self.panelCounter=0
		self.wakePanels=[]
		self.COR=np.zeros(3)
		self.sclReynolds=1.
		self.nu=0.1
		self.vortexRC=0.05
		self.rho=1.
		self.liftData=AeroSurfFactory()
		self.dragData=AeroSurfFactory()
		self.momentData=AeroSurfFactory()
	
	def getFreeVelocity(self):
		return self.VEL
	
	def setFreeVelocity(self,VEL):
		self.VEL=VEL
		
	def panelCounterIncrease(self):
		self.panelCounter+=1

	def buildGids(self):
		for r in self.getRegionList():
			for i in self.getRegion(r).getPanelList():
				self.getRegion(r).getPanel(i).connectFlowDomain(self)
				gid=self.getRegion(r).getPanel(i).gid				
				self.gids[gid]=self.getRegion(r).getPanel(i)

		
	def getPanelByGid(self,gid):
		return self.gids[gid]
				
	def setRefPT(self,PT):
		self.refPT=PT

	def addRegion(self,reg):
		self.regions[self.nReg]=reg
		self.nPan=self.nPan+self.regions[self.nReg].getNumberOfPanels()
		self.nReg+=1
		

	def getPanelGidsByBoundingBox(self,A,B,):
		l=[]
		for r in self.getRegionList():
			nPan=self.getRegion(r).getNumberOfPanels()
			ofs=self.getPanelOffset(r)
			for i in range(0,nPan):	
				p=self.getRegion(r).getPanel(i)
				c1=True
				c2=True
				for j in range(0,3):
					c1*=np.any(p.pts[:,j]>=A[j]) and np.any(p.pts[:,j]<=B[j])
					c2*=np.any(p.pts[:,j]<=A[j]) and np.any(p.pts[:,j]>=B[j])
				#if (p[0] > A[0] and p[0] < B[0]) and (p[1] > A[1] and p[1] < B[1]) and (p[2] > A[2] and p[2] < B[2]):
				if c1 or c2:
					l.append(p.gid)
		return l

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



	def reinitPanelsGeometry(self):
		for r in self.getRegionList():
			for i in self.getRegion(r).getPanelList():
				self.getRegion(r).getPanel(i).initGeometry()

	def setPanelsFreeStream(self,v):
		for r in self.getRegionList():
			for i in self.getRegion(r).getPanelList():
				self.getRegion(r).getPanel(i).setFreeStream(v)

	def getPanelOffset(self,i):
		ofs=0
		for i in range(0,i):
			ofs=ofs+self.getRegion(i).getNumberOfPanels()
		return ofs

					
	def getVelocityAtPoints(self,pts):
		pass

	def relaxWake(self,w,restricted,KW=0.05):
		pass
	
	def getForce(self):
		return self.FORCE.sum(axis=0)
	
	def getMoment(self):
		return self.MOMENT
	
	def allocMatrices(self):
		n=self.panelCounter
		self.AIB=np.zeros((n,n))
		self.AIW=np.zeros((n,n))
		self.GAMA=np.zeros(n)
		self.RHS=np.zeros(n)
		self.AREAS=np.zeros(n)
		self.CHORDS=np.zeros(n)
		self.PTS=np.zeros((n,3))
		self.VORTEXPTS=np.zeros((n*4,3))
		self.NORMALS=np.zeros((n,3))
		self.TANGENTS=np.zeros((n,3))
		self.GAMAV=np.zeros((n,3))
		self.GAMAVNorm=np.zeros(n)
		self.GAMAVUNIT=np.zeros((n,3))
		self.GAMA=np.empty(n, dtype=np.double)
		self.RHS=np.zeros(n, dtype=np.double)
		self.FORCE=np.zeros((n,3), dtype=np.double)
		self.FORCEARM=np.zeros((n,3), dtype=np.double)
		self.MOMENT=np.zeros(3)
		self.ALPHA=np.zeros(n, dtype=np.double)
		self.ALPHAEFF=np.zeros(n, dtype=np.double)
		self.CL=np.empty(n, dtype=np.double)
		self.DGAMA=np.zeros(n)
		self.GAMAGRAD=np.zeros(n)
		self.GAMARES=np.zeros(n)
		self.DECAMBER=np.zeros(n)
		self.useVisc=np.zeros(n,dtype=np.bool)
		self.useGamaGrad=np.zeros(n,dtype=np.bool)

	def allocWakePoints(self):
		self.AIW.fill(0.0)
		m=len(self.AIB)
		n=m*4##len(self.VORTEXPTS)
		#self.WWMM=[]
		nwpt=0
		nwpt_max=0
		nbp=len(self.panels)
		self.WAKEMASK=np.zeros(nbp,dtype=np.bool)
		for j in self.panels:
			ppb=self.panels[j]
			if ppb.isLast:
				self.WAKEMASK[j]=True
				nwpt+=len(ppb.wakes)
				if len(ppb.wakes) > nwpt_max:
					nwpt_max=len(ppb.wakes)
		
		self.WAKEPTS=np.zeros((nbp,nwpt_max,4,3))
		for j in self.panels:
			if self.WAKEMASK[j]:
				ppb=self.panels[j]
				wpn=0
				for w in ppb.wakes:
					pp=ppb.wakes[w]
					for k in range(0,4):
						self.WAKEPTS[j][wpn][k]=pp.pts[k]
					wpn+=1
					
						
	def computeWakePts(self):
		for i in self.regions:
			k=0
			wakes=self.regions[i].wakes
			n=len(wakes)
			nc=len(wakes[list(wakes.keys())[0]]['panels'])
			self.regions[i].wakeMask=list(wakes.keys())
			self.regions[i].WAKEPTS=np.zeros((4*n*nc,3))
			for j in self.regions[i].wakes:#['panels']:
				for jj in self.regions[i].wakes[j]['panels']:
					self.regions[i].WAKEPTS[k*4:(k+1)*(4)]=self.regions[i].wakes[j]['panels'][jj].pts
					k+=1
			self.regions[i].nWakePanels=k
			
		
	def assembleMatrices(self):
		k=0
		self.panels={}
		for i in self.regions:
			for j in self.regions[i].panels:
				self.panels[k]=self.regions[i].panels[j]
				self.PTS[k]=self.regions[i].panels[j].getCpoint()
				self.VORTEXPTS[k*4:(k+1)*(4)]=self.regions[i].panels[j].pts
				self.NORMALS[k]=self.regions[i].panels[j].normalVector
				self.TANGENTS[k]=self.regions[i].panels[j].tangentVector
				pts=self.regions[i].panels[j].pts
				gv=pts[1]-pts[0]
				self.GAMAV[k]=gv
				self.GAMAVNorm[k]=np.linalg.norm(gv)
				self.GAMAVUNIT[k]=gv/self.GAMAVNorm[k]
				self.AREAS[k]=self.regions[i].panels[j].getArea()
				self.CHORDS[k]=self.regions[i].panels[j].chord

				k+=1
	
	def assembleMatrixA_numba(self):
		assembleMatrixA_numba(self.AIB,self.PTS,self.VORTEXPTS,self.NORMALS)
		
	def computeWakeContribution_numba(self):
		computeWakeContribution_numba(self.AIW,self.PTS,self.WAKEPTS,self.WAKEMASK,self.NORMALS)
	
	def assembleMatrixA(self):
		m=len(self.AIB)
		n=m*4##len(self.VORTEXPTS)
		
		ia=np.linspace(0,n-1,n,dtype=np.int)
		ib=np.linspace(1,n,n,dtype=np.int); ib[3::4]=ia[0::4];
		
		for i in range(0,m):
			P=self.PTS[i]
			a=P-self.VORTEXPTS[ia]
			b=P-self.VORTEXPTS[ib]
			an=np.linalg.norm(a,axis=1)
			bn=np.linalg.norm(b,axis=1)
			K=(an+bn)/(an*bn*( an*bn+np.sum(np.multiply(a,b),axis=1)))/4./np.pi
			V=np.cross(a,b,axis=1) *  K[np.newaxis].T
			#self.AIB[i] =  np.sum( np.multiply( V[0::4]+V[1::4]+V[2::4]+V[3::4], self.NORMALS[i]), axis=1)
			self.AIB[i] =  np.dot( V[0::4]+V[1::4]+V[2::4]+V[3::4], self.NORMALS[i])

	
	def assembleMatrixA_slow(self):
		m=len(self.AIB)
		n=m*4##len(self.VORTEXPTS)
		self.AIB.fill(0.0)
		for i in range(0,m):
			P=self.PTS[i]
			for j in self.panels:
				pp=self.panels[j]
				aa=0.0
				for k in range(0,4):
					a=P-pp.pts[k]
					b=P-pp.pts[(k+1)%4]
					an=np.linalg.norm(a)
					bn=np.linalg.norm(b)
					#K= (A+B)/(A*B*(A*B+vect_dot(a,b) ))/4.0/M_PI;
					K=(an+bn)/(an*bn*( an*bn+np.dot(a,b) ))/4./np.pi # np.sum(np.multiply(a,b))))
					V=np.cross(a,b) *  K#[np.newaxis].T
					self.AIB[i][j]+=  np.sum( np.multiply( V, self.NORMALS[i]))
					aa+=  np.sum( np.multiply( V, self.NORMALS[i]))
					#print i,j,a,b,P,V,self.NORMALS[i], np.sum( np.multiply( V, self.NORMALS[i]))
					print(i,j,V,self.NORMALS[i], np.sum( np.multiply( V, self.NORMALS[i])),np.dot(V,self.NORMALS[i]))
				print('======================\n')
				print( aa, '\n\n')

	
	def assembleWakeContribution(self):
		nw=0
		nt=0
		for r in self.regions:
			nw+=self.regions[r].nWakePanels
			nt+=self.regions[r].nSpanwise
			
		self.wakeNORMALS=np.zeros((nw,3))
		self.wakeMASK=np.zeros(nt,dtype=np.int)
		#self.wakeMASKB=np.zeros(nw*4,dtype=bool)
		self.wakeMASKB=np.zeros(self.nPan,dtype=bool)
		self.wakeStripMASKI=np.zeros((self.nPan,nw*4),dtype=int)
		self.wakeStripMASKB=np.zeros((self.nPan,self.nPan*4),dtype=bool)

		self.wakeIA=np.zeros(nw*4,dtype=np.int)
		self.wakeIB=np.zeros(nw*4,dtype=np.int)
		self.wakePTS=np.zeros((nw*4,3))
			
		k=0
		kk=0
		kt=0
		for r in self.regions:
			nt=self.regions[r].nSpanwise
			nn=self.regions[r].nWakePanels
			n=self.regions[r].nWakePanels*4
			
			
			ia=np.linspace(k,k+n-1,n,dtype=np.int)
			ib=np.linspace(k+1,k+n,n,dtype=np.int); ib[3::4]=ia[0::4];
			self.wakeIA[k:k+n]=ia
			self.wakeIB[k:k+n]=ib
			#self.wakeNORMALS[kk:kk+nn]=self.NORMALS[self.regions[r].wakeMask]
			self.wakeMASK[kt:kt+nt]=self.regions[r].wakeMask
			self.wakeMASKB[list(np.array(self.regions[r].wakeMask,dtype=int)) ]=True
			#self.wakeMASKB[list(np.array(self.regions[r].wakeMask,dtype=int) *4+0)]=True
			#self.wakeMASKB[list(np.array(self.regions[r].wakeMask,dtype=int) *4+1)]=True
			#self.wakeMASKB[list(np.array(self.regions[r].wakeMask,dtype=int) *4+2)]=True
			#self.wakeMASKB[list(np.array(self.regions[r].wakeMask,dtype=int) *4+3)]=True
			self.wakePTS[k:k+n]=self.regions[r].WAKEPTS
			k+=n
			kk+=nn
			kt+=nt
		
		ii=0
		
		for r in self.regions:
			nt=self.regions[r].nSpanwise
			nn=self.regions[r].nWakePanels
			nch=int(nn/nt)
			n=self.regions[r].nWakePanels*4
			nwp=nch*4
			for pi in self.regions[r].panels:
				jj=0
				for i in self.regions[r].panels:
					p=self.regions[r].panels[i]
					if p.isLast:
						self.wakeStripMASKB[ii][i*4:i*4+4]=True
					jj+=nwp
				ii+=1


	def computeWakeContribution_slow(self):
		self.AIW.fill(0.0)
		m=len(self.AIB)
		n=m*4##len(self.VORTEXPTS)
		#self.WWMM=[]
		for i in range(0,m):
			P=self.PTS[i]
			for j in self.panels:
				ppb=self.panels[j]
				if ppb.isLast:
					print(i)
					#self.WWMM.append(j)
					for w in ppb.wakes:
						pp=ppb.wakes[w]
						for k in range(0,4):
							a=P-pp.pts[k]
							b=P-pp.pts[(k+1)%4]
							an=np.linalg.norm(a)
							bn=np.linalg.norm(b)
							#K= (A+B)/(A*B*(A*B+vect_dot(a,b) ))/4.0/M_PI;
							K=(an+bn)/(an*bn*( an*bn+np.dot(a,b) ))/4./np.pi # np.sum(np.multiply(a,b))))
							V=np.cross(a,b) *  K#[np.newaxis].T
							self.AIW[i][j]+=  np.sum( np.multiply( V, self.NORMALS[i]))
	
	def computeWakeContribution(self):
		self.AIW.fill(0.0)
		m=len(self.AIB)
		rc=0.001
		for i in range(0,m):
		#for i in self.wakeMASK:
			#if i in self.wakeMASK:
			print(i)
			P=self.PTS[i]
			#a=P-self.wakePTS[self.wakeIA][self.wakeMASKB]
			#b=P-self.wakePTS[self.wakeIB][self.wakeMASKB]
			#c=self.wakePTS[self.wakeIB][self.wakeMASKB]-self.wakePTS[self.wakeIA][self.wakeMASKB]
			a=P-self.wakePTS[self.wakeIA]
			b=P-self.wakePTS[self.wakeIB]
			c=self.wakePTS[self.wakeIB]-self.wakePTS[self.wakeIA]
			rd=np.linalg.norm(np.cross(a,b),axis=1)/np.linalg.norm(c)
			
			an= np.sum(a*a,axis=1)**0.5#np.linalg.norm(a,axis=1)
			bn=np.sum(b*b,axis=1)**0.5#np.linalg.norm(b,axis=1)
			K=(an+bn)/(an*bn*( an*bn+np.sum(np.multiply(a,b),axis=1)))/4./np.pi
			####K*=(1-np.exp(-(rd/rc**2)))*2.
			# BUG!!! viscous core not working
			#K[rd<rc]*=rd[rd<rc]*rd[rd<rc]/rc/rc
			V=np.cross(a,b,axis=1) *  K[np.newaxis].T
			xx=np.sum( np.multiply( V[0::4]+V[1::4]+V[2::4]+V[3::4],self.NORMALS[i]), axis=1)
			x=np.split(xx,len(self.wakeMASK))
			for j in range(0,len(x)):
				self.AIW[i][self.wakeMASK[j]]+=np.sum(x[j])
			#self.AIW[i][self.wakeMASK]+=np.sum( np.multiply( V[0::4]+V[1::4]+V[2::4]+V[3::4],self.NORMALS[i]), axis=1)
			#self.AIW[i][self.wakeMASK]+=np.dot( V[0::4]+V[1::4]+V[2::4]+V[3::4],self.wakeNORMALS[i])

	def computeWakeContributionBIG(self):
		#self.AIW.fill(0.0)
		a=self.big_PTS-self.big_wakePTS[self.big_wakeIA]
		b=self.big_PTS-self.big_wakePTS[self.big_wakeIB]
		an= np.sum(a*a,axis=1)**0.5#np.linalg.norm(a,axis=1)
		bn=np.sum(b*b,axis=1)**0.5#np.linalg.norm(b,axis=1)
		K=(an+bn)/(an*bn*( an*bn+np.sum(np.multiply(a,b),axis=1)))/4./np.pi
		#V=np.cross(a,b,axis=1) *  K[np.newaxis].T
		
		self._vv[:,0]=a[:,1]*b[:,2] -a[:,2]*b[:,1]
		self._vv[:,1]=(a[:,0]*b[:,2] -a[:,2]*b[:,0])*-1
		self._vv[:,2]=a[:,0]*b[:,1] -a[:,1]*b[:,0]
		V= K.reshape(self.big_n,1) * self._vv
		AIWAKE =  np.sum( np.multiply( V[0::4]+V[1::4]+V[2::4]+V[3::4],self.big_wakeNORMALS), axis=1)
		#self.AIW[i][self.wakeMASK]+=AIWAKE
		self.AIW=AIWAKE.reshape(self.nPan,self.nPan).T
	
	#def computeWakeContributionCY(self):
		#def _computeWakeContributionCY(np.ndarray[np.float64_t, ndim=2] AIW,mm,nn,np.ndarray[np.float64_t, ndim=2] PTS,
				 #np.ndarray[np.float64_t, ndim=2] wakePTS,
				 #np.ndarray[np.int_t, ndim=1] wakeIA,
				 #np.ndarray[np.int_t, ndim=1] wakeIB,
				 #np.ndarray[np.float64_t, ndim=2] wakeMASK,
				 #np.ndarray[np.float64_t, ndim=2] wakeNORMALS):
			
			#cdef Py_ssize_t i
			#cdef int m=mm
			#cdef int n=nn
			
			##np.ndarray[np.float64_t, ndim=1]
	
			#cdef np.ndarray[np.float64_t, ndim=2]vv=np.zeros((n,3))
			#for i in xrange(0,m):
				#P=PTS[i]
				#a=P-wakePTS[wakeIA]
				#b=P-wakePTS[wakeIB]
				#an= np.sum(a*a,axis=1)**0.5
				#bn=np.sum(b*b,axis=1)**0.5
				#K=(an+bn)/(an*bn*( an*bn+np.sum(np.multiply(a,b),axis=1)))
				##V=np.cross(a,b,axis=1) *  K[np.newaxis].T
				#vv[:,0]=a[:,1]*b[:,2] -a[:,2]*b[:,1]
				#vv[:,1]=(a[:,0]*b[:,2] -a[:,2]*b[:,0])*-1
				#vv[:,2]=a[:,0]*b[:,1] -a[:,1]*b[:,0]
				
				#V= K.reshape(n,1) * vv
				#AIWAKE =  np.sum( np.multiply( V[0::4]+V[1::4]+V[2::4]+V[3::4],wakeNORMALS), axis=1)
				#AIW[i][wakeMASK]+=AIWAKE
				
		#self.AIW.fill(0.0)
		#m=len(self.AIW)
		#n=len(self.wakeIA)
		#_computeWakeContributionCY(self.AIW,m,n,self.PTS,self.wakePTS,self.wakeIA,self.wakeIB,self.wakeMASK,self.wakeNORMALS)
	
	#def computeWakeContributionNB(self):
		##@nb.jit(nopython=True)	
		#def _computeWakeContributionNB(AIW,m,n,PTS,wakePTS,wakeIA,wakeIB,wakeMASK,wakeNORMALS):
			#vv=np.zeros((n,3))
			#for i in xrange(0,m):
				#P=PTS[i]
				#a=P-wakePTS[wakeIA]
				#b=P-wakePTS[wakeIB]
				#an= np.sum(a*a,axis=1)**0.5
				#bn=np.sum(b*b,axis=1)**0.5
				#K=(an+bn)/(an*bn*( an*bn+np.sum(np.multiply(a,b),axis=1)))
				##V=np.cross(a,b,axis=1) *  K[np.newaxis].T
				#vv[:,0]=a[:,1]*b[:,2] -a[:,2]*b[:,1]
				#vv[:,1]=(a[:,0]*b[:,2] -a[:,2]*b[:,0])*-1
				#vv[:,2]=a[:,0]*b[:,1] -a[:,1]*b[:,0]
				
				#V= K.reshape(n,1) * vv
				#AIWAKE =  np.sum( np.multiply( V[0::4]+V[1::4]+V[2::4]+V[3::4],wakeNORMALS), axis=1)
				#AIW[i][wakeMASK]+=AIWAKE
				
		#self.AIW.fill(0.0)
		#m=len(self.AIW)
		#n=len(self.wakeIA)
		#_computeWakeContributionNB(self.AIW,m,n,self.PTS,self.wakePTS,self.wakeIA,self.wakeIB,self.wakeMASK,self.wakeNORMALS)
	
	
	
	def initSolution(self):
		self.allocMatrices()
		self.assembleMatrices()
		self.allocWakePoints()
		self.assembleMatrixA_numba()
		self.computeWakePts()
		self.assembleWakeContribution()
		self.computeWakeContribution_numba()
		self.FORCEARM=self.PTS-self.COR
		self.NORMALS_X=self.NORMALS.copy()

	def updateTimeStep(self):
		#self.updateGeometry()
		if self.fcs is not None:
			sefl.fcs.update()
		self.assembleMatrices()
		#self.allocWakePoints()
		self.assembleMatrixA_numba()
		self.computeWakePts()
		self.assembleWakeContribution()
		self.computeWakeContribution_numba()
		self.FORCEARM=self.PTS-self.COR
		#self.NORMALS_X=self.NORMALS.copy()
	
	
	def updateGeometry(self):
		#for i in self.regions:
			#wakes=self.regions[i].wakes
			#n=len(wakes)
			#nc=len(wakes[wakes.keys()[0]]['panels'])
			#self.regions[i].wakeMask=wakes.keys()
			#for j in self.regions[i].wakes:
				#for jj in self.regions[i].wakes[j]['panels']:
					#self.modifyPanelGeometry(self.regions[i].wakes[j]['panels'][jj])
			
			#for j in self.regions[i].panels:
				#self.modifyPanelGeometry(self.regions[i].panels[j])			
				
		k=0
		self.panels={}
		for i in self.regions:
			for j in self.regions[i].panels:
				#self.panels[k]=self.regions[i].panels[j]
				#p=self.regions[i].panels[j]
				self.regions[i].panels[j].controlPoint=self.PTS[k]
				self.regions[i].panels[j].pts=self.VORTEXPTS[k*4:(k+1)*(4)]
				k+=1		
		self.reinitPanelsGeometry()
		
		
	def computeVelocitiesTotal(self):
		self.vref=np.linalg.norm(self.VEL)
		self.VTOT=np.tile(self.VEL,(self.nPan,1)) + np.cross(self.OMEGA,self.PTS-self.COR)
		self.vnorm=np.linalg.norm(self.VTOT,axis=1)
	
	#@jit(nopython=True,nogil=True,parallel=True)
	def compute(self):
		self.computeVelocitiesTotal()
		self.RHS=np.sum(-self.VTOT* self.NORMALS,axis=1)
		self.GAMA_OLD=self.GAMA.copy()
		self.GAMA=np.linalg.solve(self.AIB+self.AIW,self.RHS)
		if np.all(self.GAMA_OLD == 0):
			self.dGAMA_dt=np.zeros_like(self.GAMA)
		else:
			if hasattr(self,'dt'):
				self.dGAMA_dt=(self.GAMA-self.GAMA_OLD)/self.dt
			else:
				self.dGAMA_dt=np.zeros_like(self.GAMA)
				
		#self.ALPHA=self.GAMA/(self.vnorm*0.5 *2.*np.pi)#*self.AREAS)
		
		#self.dLift=self.GAMA.copy()
		vref=np.linalg.norm(self.vref)
		
		k=0
		self.FORCE=np.zeros((self.nPan,3))# np.cross(self.VTOT,self.GAMAV)
		for r in self.regions:
			region=self.regions[r]
			m=region.nSpanwise 
			n=region.nChordwise
			for i in range(0,m):
				p=region.panels[i*n]
				p.G=self.GAMA[k]
				pn=region.panels[i*n+1]
				pn.G=self.GAMA[k+1]
				#p.dp=np.dot(self.VTOT[k],p.tangentVector_unit)*(p.G-pn.G)/(p.chord/2 +pn.chord/2)#+self.dGAMA_dt[k]
				#self.FORCE[k]=p.normalVector*p.dp*p.getArea()
				self.FORCE[k]=p.normalVector* p.G*self.GAMAVNorm[k] *self.rho *np.dot(p.tangentVector_unit,self.VTOT[k])# + self.dGAMA_dt[k]
				if not np.all(np.isfinite(self.FORCE[k])):
					self.FORCE[k].fill(0)
				p.forceVector=self.FORCE[k]
				sign=np.sign(np.dot(p.forceVector,p.normalVector))
				p.pressureCoef=sign*np.linalg.norm(p.forceVector)/vref/p.getArea()/self.rho
				k+=1
				for j in range(1,n):
					p=region.panels[i*n+j]
					p.G=self.GAMA[k]
					pp=region.panels[i*n+j-1]
					#p.dp=np.dot(self.VTOT[k],p.tangentVector_unit)*(p.G-pp.G)/(p.chord/2 +pp.chord/2)#+self.dGAMA_dt[k]
					#self.FORCE[k]=p.normalVector*p.dp*p.getArea()
					self.FORCE[k]=p.normalVector*self.GAMAVNorm[k]* (p.G-pp.G)*self.rho *np.dot(p.tangentVector_unit,self.VTOT[k])# + self.dGAMA_dt[k]
					if not np.all(np.isfinite(self.FORCE[k])):
						self.FORCE[k].fill(0)
					p.forceVector=self.FORCE[k]
					p.pressureCoef=np.linalg.norm(p.forceVector)/vref/p.getArea()/self.rho
					k+=1
		self.MOMENT=np.sum(np.cross(self.FORCEARM,self.FORCE),axis=0)
		#k=0
		#for i in self.regions:
			#for j in self.regions[i].panels:
				#self.dLift[k]=self.regions[i].panels[j].dLift
				#k+=1
		#self.FORCE=np.cross(self.VEL,self.GAMAV)*self.dLift.reshape(self.nPan,1)
		#self.CL=self.dLift/self.AREAS/self.vref**2.
		
	def computeViscous(self,maxIt=1):
		self.computeVelocitiesTotal()
		#self.DECAMBER.fill(0)
		self.CL_OLD=np.zeros_like(self.CL)
		self.ALPHACOR=np.zeros_like(self.ALPHA)
		for i in range(0,maxIt):
			for j in range(0,self.nPan):
				self.NORMALS[j]=bc.rotAx(np.array([self.NORMALS_X[j]]),self.GAMAVUNIT[j],np.zeros(3),self.DECAMBER[j])
			self.assembleMatrixA()
			
			self.RHS= np.sum(-self.VTOT* self.NORMALS,axis=1)
			self.GAMA=np.linalg.solve(    self.AIB+self.AIW , self.RHS )

			self.CL=self.GAMA*self.GAMAVNorm/(self.vnorm *0.5 *self.AREAS)
			#self.ALPHA=self.GAMA/(self.vnorm*0.5 *2.*np.pi*self.AREAS)	
			self.ALPHA=self.CL/2/np.pi
			self.ALPHAEFF=self.ALPHA-self.DECAMBER #+ (self.CL_OLD-self.CL)/2./np.pi
			#self.ALPHAEFF=self.ALPHA-self.ALPHACOR
			self.viscDataCompute()
			

			self.RESIDUAL=self.CL -self.viscCL - self.viscCL_correction
			#self.ALPHACOR=(self.CL-self.CL_OLD)/2/np.pi +self.DECAMBER
			
			self.DECAMBER=self.RESIDUAL/2./np.pi
		
			self.CL_OLD=self.CL.copy()
		
		self.GAMAGRAD=np.ediff1d(self.GAMA,to_begin=0)
		self.GAMARES=self.GAMA.copy()
		self.GAMARES[self.useGamaGrad]=self.GAMAGRAD[self.useGamaGrad]

		#self.DGAMA=self.GAMA/self.GAMAVNorm*(self.vnorm *0.5 *self.AREAS)
		#self.GAMARES[self.useVisc]+=self.DGAMA[self.useVisc]

		forceLift=0.5*self.vnorm*self.vnorm*self.AREAS*self.CL
		forceDrag=0.5*self.vnorm*self.vnorm*self.AREAS*self.viscCD
		vnorm=self.VTOT/np.linalg.norm(self.VTOT,axis=1).reshape(self.nPan,1)
		self.FORCE=np.cross(vnorm,self.GAMAVUNIT)*forceLift.reshape(self.nPan,1)
		self.FORCE_DRAG=vnorm*forceDrag.reshape(self.nPan,1)
		self.FORCE+=self.FORCE_DRAG
		mmag=0.5*self.vnorm*self.vnorm*self.AREAS*self.viscCM*self.CHORDS
		mmm=self.GAMAVUNIT*mmag.reshape(self.nPan,1)
		self.MOMENT=np.sum(np.cross(self.FORCEARM,self.FORCE),axis=0) + np.sum(mmm[self.useVisc],axis=0)
	
	

	def computeViscousStep(self):
		self.computeVelocitiesTotal()
		for j in range(0,self.nPan):
			self.NORMALS[j]=bc.rotAx(np.array([self.NORMALS_X[j]]),self.GAMAVUNIT[j],np.zeros(3),self.DECAMBER[j])
		self.assembleMatrixA()
		
		self.RHS= np.sum(-self.VTOT* self.NORMALS,axis=1)
		self.GAMA=np.linalg.solve(    self.AIB+self.AIW , self.RHS )

		self.CL=self.GAMA*self.GAMAVNorm/(self.vnorm *0.5 *self.AREAS)
		self.ALPHA=self.CL/2/np.pi+self.DECAMBER#(self.GAMA)/(self.vnorm*0.5 *2.*np.pi)#*self.AREAS)
		self.viscDataCompute()
		
		self.KK=2*np.pi*np.ones_like(self.CL) #self.CL/self.ALPHA
		#self.KK=(self.viscCL + self.viscCL_correction)/self.ALPHA

		self.RESIDUAL=self.viscCL + self.viscCL_correction - self.CL
		##self.DECAMBER=self.RESIDUAL/self.KK

		self.GAMAGRAD=np.ediff1d(self.GAMA,to_begin=0)
		self.GAMARES=self.GAMA.copy()
		self.GAMARES[self.useGamaGrad]=self.GAMAGRAD[self.useGamaGrad]

		forceLift=0.5*self.vnorm*self.vnorm*self.AREAS*self.CL
		forceDrag=0.5*self.vnorm*self.vnorm*self.AREAS*self.viscCD
		vnorm=self.VTOT/np.linalg.norm(self.VTOT,axis=1).reshape(self.nPan,1)
		self.FORCE=np.cross(vnorm,self.GAMAVUNIT)*forceLift.reshape(self.nPan,1)
		self.FORCE_DRAG=vnorm*forceDrag.reshape(self.nPan,1)
		self.FORCE+=self.FORCE_DRAG
		mmag=0.5*self.vnorm*self.vnorm*self.AREAS*self.viscCM*self.CHORDS
		mmm=self.GAMAVUNIT*mmag.reshape(self.nPan,1)
		self.MOMENT=np.sum(np.cross(self.FORCEARM,self.FORCE),axis=0) + np.sum(mmm[self.useVisc],axis=0)



	def computeViscousResidual(self,X):
		#X = dAlpha
		self.computeVelocitiesTotal()
		for j in range(0,self.nPan):
			self.NORMALS[j]=bc.rotAx(np.array([self.NORMALS_X[j]]),self.GAMAVUNIT[j],np.zeros(3),X[j])
			#pass
		self.assembleMatrixA()
		
		self.RHS= np.sum(-self.VTOT* self.NORMALS,axis=1)
		#self.GAMA=np.linalg.solve(    (self.AIB+self.AIW)*np.cos(X).reshape(self.nPan,1) , self.RHS*np.cos(X) )
		self.GAMA=np.linalg.solve(    self.AIB+self.AIW , self.RHS )
		self.CL=self.GAMA*self.GAMAVNorm/(self.vnorm *0.5 *self.AREAS)
		self.ALPHA=self.CL/2/np.pi
		self.ALPHAEFF=self.ALPHA-X
		self.viscDataCompute()
		

		RESIDUAL=self.CL - self.viscCL - self.viscCL_correction
		return RESIDUAL
	
	def computeViscous1(self,xtol=1e-9):
		def findResidualCL(x):
			return self.computeViscousResidual(x)
		self.DECAMBER.fill(0.)
		self.computeViscous()	
		
		sol=optimize.fsolve(findResidualCL ,self.DECAMBER,xtol=xtol)#,method='krylov')
		#sol=optimize.root(findResidualCL ,self.DECAMBER)
		self.DECAMBER=sol
		
		self.RESIDUAL=findResidualCL(self.DECAMBER)
		
		#self.computeViscousStep()
		self.GAMAGRAD=np.ediff1d(self.GAMA,to_begin=0)
		self.GAMARES=self.GAMA.copy()
		self.GAMARES[self.useGamaGrad]=self.GAMAGRAD[self.useGamaGrad]

		forceLift=0.5*self.vnorm*self.vnorm*self.AREAS*self.viscCL
		forceDrag=0.5*self.vnorm*self.vnorm*self.AREAS*self.viscCD
		vnorm=self.VTOT/np.linalg.norm(self.VTOT,axis=1).reshape(self.nPan,1)
		self.FORCE=np.cross(vnorm,self.GAMAVUNIT)*forceLift.reshape(self.nPan,1)
		#self.FORCE=self.NORMALS*forceLift.reshape(self.nPan,1)

		self.FORCE_DRAG=vnorm*forceDrag.reshape(self.nPan,1)
		self.FORCE+=self.FORCE_DRAG
		
		self.W_IND=np.dot(self.AIW,self.GAMA)
		self.FORCE_DRAG_IND=self.GAMA*self.GAMAVNorm*self.W_IND*0.25## PICOVINA!!!!!!!!!! 
		
		self.FORCE+=vnorm*self.FORCE_DRAG_IND.reshape(self.nPan,1)
		
		mmag=0.5*self.vnorm*self.vnorm*self.AREAS*self.viscCM*self.CHORDS
		mmm=self.GAMAVUNIT*mmag.reshape(self.nPan,1)
		self.MOMENT=np.sum(np.cross(self.FORCEARM,self.FORCE),axis=0) + np.sum(mmm[self.useVisc],axis=0)
		for j in self.panels:
			self.panels[j].G=self.GAMA[j]
			self.panels[j].pressureCoef=forceLift[j]# self.viscCL[j]
		
	def viscDataInit(self,alphaRng=None):
		self.viscInterp=np.linspace(0,self.nPan-1,self.nPan)
		self.liftData.setYi(self.viscInterp)
		self.dragData.setYi(self.viscInterp)
		self.momentData.setYi(self.viscInterp)
		if alphaRng is not None:
			self.liftData.alphaRng=alphaRng
			self.dragData.alphaRng=alphaRng
			self.momentData.alphaRng=alphaRng
		self.liftData.xi=np.zeros((self.nPan,2));self.liftData.xi[:,0]=self.viscInterp;
		self.dragData.xi=np.zeros((self.nPan,2));self.dragData.xi[:,0]=self.viscInterp;
		self.momentData.xi=np.zeros((self.nPan,2));self.momentData.xi[:,0]=self.viscInterp;
		self.viscCL=np.zeros(self.nPan)
		self.viscCD=np.zeros(self.nPan)
		self.viscCM=np.zeros(self.nPan)
		self.viscCL_correction=np.zeros(self.nPan)
		self.viscCD_correction=np.zeros(self.nPan)
		self.viscCM_correction=np.zeros(self.nPan)
		
	def viscDataMake(self):
		self.liftData.make()
		self.dragData.make()
		self.momentData.make()
	
	def viscDataComputeCL(self):
		self.liftData.xi[:,1]=self.ALPHAEFF;self.dragData.xi[:,1]=self.ALPHAEFF;self.momentData.xi[:,1]=self.ALPHAEFF;
		self.viscCL=self.liftData.getData(self.liftData.xi)#+self.viscCL_correction

	def viscDataCompute(self):
		self.viscDataComputeCL()
		self.viscCD=self.dragData.getData(self.dragData.xi)+self.viscCD_correction
		self.viscCM=self.momentData.getData(self.momentData.xi)+self.viscCM_correction
			
			
class VlmProblem(object):
	
	def __init__(self):
		try:
			del(self.panelGroups)
		except:
			pass
		try:
			del(self.dom1)
		except:
			pass
		
		self.vel=np.zeros(3)
		self.velMagn=1.0
		self.alpha=0.0
		self.beta=0.0
		self.pressure=1.0	
		self.dom1=VLMDomain()
		self.panelGroups={}

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
	
	def addRegion(self,reg,name=None):
		self.dom1.addRegion(reg)
		if name is not None:
			self.panelGroups[name]=reg
		else:
			self.panelGroups[len(self.panelGroups)]=reg
		
	def getSolution(self):
		self.dom1.setPanelsFreeStream(self.getFreeVelocity())
		self.dom1.reinitPanelsGeometry()
		#self.dom1.setPanelsG(1.0)
		self.dom1.solveAij()


			
	def plotResults(self,redraw=True,lineScaling=1.0,showAxisSystems=False,plotWake=False):
		if redraw:
			self.view=PanelViewer3D()
			if showAxisSystems:
				self.view.showAxisSystems=True
				self.view.drawAxisSystems()

		dom1=self.dom1
		for r in dom1.getRegionList():
			nPan=dom1.getRegion(r).getNumberOfPanels()
			lc=self.view.newPolyLine()
			for i in range(0,nPan):
				panelI=dom1.getRegion(r).getPanel(i)
				pp=panelI.getCpoint()
				
				if showAxisSystems:
					self.view.addCSYS(pp,panelI.csys,scaling=0.05)

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
		if showAxisSystems:
			self.view.drawAxisSystems()
		#self.view.drawLines()
		self.view.drawQuads()
#		self.view.drawPolyLines()
		self.view.drawText()
		self.view.addWidgets({'showAxes':False,'showScalarBar':True,'showXYPlane':False,'showYZPlane':False,'showZXPlane':False})
		self.view.show()


	def plotToFile(self,fname,plotWake=True):
		self.view=PanelViewer3D()
		dom1=self.dom1
		for r in dom1.getRegionList():
			reg=dom1.getRegion(r)
			nPan=reg.getNumberOfPanels()
			for i in range(0,nPan):
				panelI=reg.getPanel(i)
				#if not panelI.isWake or plotWake:
				self.view.addVlmPanel(panelI)
			if plotWake:
				for iw in reg.wakes:
					for jw in reg.wakes[iw]['panels']:
						self.view.addVlmPanel(reg.wakes[iw]['panels'][jw])
					
		self.view.drawQuads()
		self.view.writeVTK(fname)

	def saveAeroToFile(self,fname):
		n=len(self.dom1.gids)
		d=np.zeros((n,8))
		ks=self.dom1.gids.keys()
		for i in range(0,n):
			p=self.dom1.getPanelByGid(ks[i])
			d[i,0]=ks[i]
			d[i,1]=p.S
			d[i,2]=p.G
			d[i,3]=p.getPressureCoef()
			d[i,4]=p.dP
			d[i,5:8]=p.forceVector
		
		np.savetxt(fname,d)
		
	def loadAeroFromFile(self,fname):
		d=np.loadtxt(fname)
		nn=d.shape[0]
		n=self.dom1.nPan
		if n != nn:
			return
		ks=self.dom1.gids.keys()
		for j in range(0,n):
			i=int(d[j,0])
			p=self.dom1.getPanelByGid(i)
			p.S=d[j,1]
			p.G=d[j,2]
			p.pressureCoef=d[j,3]
			p.dP=d[j,4]
			p.forceVector=d[j,5:8]
	
		for r in self.dom1.getRegionList():
			nPan=self.dom1.getRegion(r).getNumberOfPanels()
			region=self.dom1.getRegion(r)
			for i in range(0,nPan):	
				region.getPanel(i).setForce()
		
	def getGidHT(self,addWakes=False):
		gidHT={}
		gids=self.dom1.gids.keys()
		j=0
		for i in range(0,len(gids)):
			p=self.dom1.getPanelByGid(gids[i])
			if not p.isWake or addWakes:
				gidHT[j]=gids[i]
				j+=1
		return gidHT
	
	def getBodyPanels(self):
		gidHT=self.getGidHT(addWakes=False)
		p={}
		for i in gidHT:
			p[i]=self.dom1.getPanelByGid(gidHT[i])
		return p
	
	
		
