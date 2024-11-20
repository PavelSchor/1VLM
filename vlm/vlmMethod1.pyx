from vlm.vlmPanels import VlmPanel
import numpy as np
cimport numpy as np


from collections import OrderedDict
import copy
from utils.transform import EulerMatrix
from utils.b3Vect import B3vect 
vect=B3vect()
from vlm.viewer import *

from time import gmtime, strftime
from scipy.sparse.linalg import lgmres

import affinity
import multiprocessing


#def sig_handler(signum, frame):
	#print('Segmentation Fault')

#signal.signal(signal.SIGSEGV, sig_handler) 

affinity.set_process_affinity_mask(0,2**multiprocessing.cpu_count()-1)


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
				r=panelI.midPoint -PT
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


class PanelGroup(SubSection):
	def addSubsection(self,r,join=True):
		for i in r.getPanelList():
			self.addPanel(r.getPanel(i))
		if join:
			self.nSpanwise+=r.nSpanwise
			self.nLPanChordwise=r.nLPanChordwise
		
	def addPanelGroup(self,pg):
		self.addSubsection(pg)

class GlobalDomain(object):
	def __init__(self):
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
		self.omega=np.zeros(3)
		self.sclReynolds=1.
		self.nu=0.1
		self.vortexRC=0.05
		
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

	def allocMatrices(self):
		n=self.panelCounter
		self.AIB=np.zeros((n,n))
		self.AIW=np.zeros((n,n))
		self.GAMA=np.zeros(n)
		self.RHS=np.zeros(n)
		self.PTS=np.zeros((n,3))
		self.VORTEXPTS=np.zeros((n*4,3))
		self.NORMALS=np.zeros((n,3))
	
	def computeWakePts(self):
		for i in self.regions:
			k=0
			wakes=self.regions[i].wakes
			n=len(wakes)
			nc=len(wakes[wakes.keys()[0]]['panels'])
			self.regions[i].wakeMask=wakes.keys()
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
				k+=1
	
	def assembleMatrixA(self):
		m=np.alen(self.AIB)
		n=m*4##np.alen(self.VORTEXPTS)
		
		ia=np.linspace(0,n-1,n,dtype=np.int)
		ib=np.linspace(1,n,n,dtype=np.int); ib[3::4]=ia[0::4];
		
		for i in range(0,m):
			P=self.PTS[i]
			a=P-self.VORTEXPTS[ia]
			b=P-self.VORTEXPTS[ib]
			an=np.linalg.norm(a,axis=1)
			bn=np.linalg.norm(b,axis=1)
			K=(an+bn)/(an*bn*( an*bn+np.sum(np.multiply(a,b),axis=1)))
			V=np.cross(a,b,axis=1) *  K[np.newaxis].T
			self.AIB[i] =  np.sum( np.multiply( V[0::4]+V[1::4]+V[2::4]+V[3::4], self.NORMALS), axis=1)
		
	def assembleWakeContribution(self):
		nw=0
		
		for r in self.regions:
			nw+=self.regions[r].nWakePanels
		
		self.wakeNORMALS=np.zeros((nw,3))
		self.wakeMASK=np.zeros(nw,dtype=np.int)
		self.wakeIA=np.zeros(nw*4,dtype=np.int)
		self.wakeIB=np.zeros(nw*4,dtype=np.int)
		self.wakePTS=np.zeros((nw*4,3))
		
		
		k=0
		kk=0
		for r in self.regions:
			nn=self.regions[r].nWakePanels
			n=self.regions[r].nWakePanels*4
			
			self.wakeNORMALS[kk:kk+nn]=self.NORMALS[self.regions[r].wakeMask]
			self.wakeMASK[kk:kk+nn]=self.regions[r].wakeMask
			ia=np.linspace(0,n-1,n,dtype=np.int)
			ib=np.linspace(1,n,n,dtype=np.int); ib[3::4]=ia[0::4];
			self.wakeIA[k:k+n]=ia
			self.wakeIB[k:k+n]=ib
			
			self.wakePTS[k:k+n]=self.regions[r].WAKEPTS
			
			k+=n
			kk+=nn
	def computeWakeContribution(self):
		self.AIW.fill(0.0)
		m=np.alen(self.AIB)
		for i in range(0,m):
		#it = np.nditer(a, flags=['f_index'])
		#while not it.finished:
			#print "%d <%d>" % (it[0], it.index),
			#it.iternext()
			P=self.PTS[i]
			a=P-self.wakePTS[self.wakeIA]
			b=P-self.wakePTS[self.wakeIB]
			an= np.sum(a*a,axis=1)**0.5#np.linalg.norm(a,axis=1)
			bn=np.sum(b*b,axis=1)**0.5#np.linalg.norm(b,axis=1)
			K=(an+bn)/(an*bn*( an*bn+np.sum(np.multiply(a,b),axis=1)))
			V=np.cross(a,b,axis=1) *  K[np.newaxis].T
			AIWAKE =  np.sum( np.multiply( V[0::4]+V[1::4]+V[2::4]+V[3::4],self.wakeNORMALS), axis=1)
			self.AIW[i][self.wakeMASK]+=AIWAKE
	
	
	def computeWakeContributionCY(self):
		def _computeWakeContributionCY(np.ndarray[np.float64_t, ndim=2] AIW,
				 np.int_t mm,
				 np.int_t nn,
				 np.ndarray[np.float64_t, ndim=2] PTS,
				 np.ndarray[np.float64_t, ndim=2] wakePTS,
				 np.ndarray[np.int_t, ndim=1] wakeIA,
				 np.ndarray[np.int_t, ndim=1] wakeIB,
				 np.ndarray[np.int_t, ndim=1] wakeMASK,
				 np.ndarray[np.float64_t, ndim=2] wakeNORMALS):
			
			cdef Py_ssize_t i
			cdef int m=mm
			cdef int n=nn
			
			#np.ndarray[np.float64_t, ndim=1]
	
			cdef np.ndarray[np.float64_t, ndim=2]vv=np.zeros((n,3))
			for i in xrange(0,m):
				P=PTS[i]
				a=P-wakePTS[wakeIA]
				b=P-wakePTS[wakeIB]
				an= np.sum(a*a,axis=1)**0.5
				bn=np.sum(b*b,axis=1)**0.5
				K=(an+bn)/(an*bn*( an*bn+np.sum(np.multiply(a,b),axis=1)))
				#V=np.cross(a,b,axis=1) *  K[np.newaxis].T
				vv[:,0]=a[:,1]*b[:,2] -a[:,2]*b[:,1]
				vv[:,1]=(a[:,0]*b[:,2] -a[:,2]*b[:,0])*-1
				vv[:,2]=a[:,0]*b[:,1] -a[:,1]*b[:,0]
				
				V= K.reshape(n,1) * vv
				AIWAKE =  np.sum( np.multiply( V[0::4]+V[1::4]+V[2::4]+V[3::4],wakeNORMALS), axis=1)
				AIW[i][wakeMASK]+=AIWAKE
				
		self.AIW.fill(0.0)
		m=np.alen(self.AIW)
		n=np.alen(self.wakeIA)
		_computeWakeContributionCY(self.AIW,m,n,self.PTS,self.wakePTS,self.wakeIA,self.wakeIB,self.wakeMASK,self.wakeNORMALS)
	
	def computeWakeContributionNB(self):
		def _computeWakeContributionNB(AIW,m,n,PTS,wakePTS,wakeIA,wakeIB,wakeMASK,wakeNORMALS):
			vv=np.zeros((n,3))
			for i in xrange(0,m):
				P=PTS[i]
				a=P-wakePTS[wakeIA]
				b=P-wakePTS[wakeIB]
				an= np.sum(a*a,axis=1)**0.5
				bn=np.sum(b*b,axis=1)**0.5
				K=(an+bn)/(an*bn*( an*bn+np.sum(np.multiply(a,b),axis=1)))
				#V=np.cross(a,b,axis=1) *  K[np.newaxis].T
				vv[:,0]=a[:,1]*b[:,2] -a[:,2]*b[:,1]
				vv[:,1]=(a[:,0]*b[:,2] -a[:,2]*b[:,0])*-1
				vv[:,2]=a[:,0]*b[:,1] -a[:,1]*b[:,0]
				
				V= K.reshape(n,1) * vv
				AIWAKE =  np.sum( np.multiply( V[0::4]+V[1::4]+V[2::4]+V[3::4],wakeNORMALS), axis=1)
				AIW[i][wakeMASK]+=AIWAKE
				
		self.AIW.fill(0.0)
		m=np.alen(self.AIW)
		n=np.alen(self.wakeIA)
		_computeWakeContributionNB(self.AIW,m,n,self.PTS,self.wakePTS,self.wakeIA,self.wakeIB,self.wakeMASK,self.wakeNORMALS)
			
			
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
		self.dom1=GlobalDomain()

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
	
	def addRegion(self,reg):
		self.dom1.addRegion(reg)
		
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
			nPan=dom1.getRegion(r).getNumberOfPanels()
			for i in range(0,nPan):
				panelI=dom1.getRegion(r).getPanel(i)
				if not panelI.isWake or plotWake:
					self.view.addVlmPanel(panelI)
		self.view.drawQuads()
		self.view.writeVTP(fname)

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
	
	
		