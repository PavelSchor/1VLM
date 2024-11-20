from panelMethod.doubletSourcePanels import doubletSourcePanel
import numpy as np
import copy
from panelMethod.env import GlobalEnvironment
from utils.transform import EulerMatrix
from utils.b3Vect import B3vect 
vect=B3vect()
from panelMethod.viewer import *
from panelMethod.boundaryLayer2D import *

from time import gmtime, strftime
from scipy.sparse.linalg import lgmres

import ctypes
from numpy.ctypeslib import ndpointer

from Cython import nogil
import affinity
import multiprocessing
import pylab as plt
from pyIBLM import IblmStrip
from scipy.interpolate import UnivariateSpline
from scipy.spatial import ConvexHull
import signal

def sig_handler(signum, frame):
	print('Segmentation Fault')

signal.signal(signal.SIGSEGV, sig_handler) 

affinity.set_process_affinity_mask(0,2**multiprocessing.cpu_count()-1)

try:
	testlib = ctypes.CDLL('panelMethod/testlib_w32.so')
except:
	pass
try:
	testlib = ctypes.CDLL('panelMethod/testlib_w64.so')
except:
	pass


try:
	testlib = ctypes.CDLL('panelMethod/testlib_lnx32.so')
except:
	pass

try:
	testlib = ctypes.CDLL('panelMethod/testlib_lnx64.so')
except:
	pass

#void getUVW_at_points(double*UVW, double*poi, int n, double*pts,int npt, int ptdim_max,int*ptdim, double*ORIGS, double*csysesLoc,double*D,double*S,double*RC){
testlib.getUVW_at_points.restype = ctypes.c_void_p;
testlib.getUVW_at_points.argtypes = [ndpointer(ctypes.c_double),
ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
ctypes.c_int32,
ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
ctypes.c_int32,
ctypes.c_int32,
ndpointer(ctypes.c_int32),
ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
ndpointer(ctypes.c_double),
ndpointer(ctypes.c_double),
ndpointer(ctypes.c_double),
ndpointer(ctypes.c_int32,ndim=2, flags='C_CONTIGUOUS')];


#void getPHI_at_points(double*PHI, double*poi, int n, double*pts,int npt, int ptdim_max,int*ptdim, double*ORIGS, double*csysesLoc,double*D,double*S){
testlib.getPHI_at_points.restype = ctypes.c_void_p;
testlib.getPHI_at_points.argtypes = [ndpointer(ctypes.c_double),
ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
ctypes.c_int32,
ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
ctypes.c_int32,
ctypes.c_int32,
ndpointer(ctypes.c_int32),
ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
ndpointer(ctypes.c_double),
ndpointer(ctypes.c_double)];

#void getAij_matrix0(double*AIJ, double*poi, int n, double*pts,int npt, int ptdim_max,int*ptdim, double*ORIGS, double*csysesLoc,double*D,double*S,int*gid,int*isFirst,int*isLast,int*isWake){

testlib.getAij_matrix0.restype = ctypes.c_void_p;
testlib.getAij_matrix0.argtypes = [ndpointer(ctypes.c_double),
ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
ctypes.c_int32,
ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
ctypes.c_int32,
ctypes.c_int32,
ndpointer(ctypes.c_int32),
ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
ndpointer(ctypes.c_double),
ndpointer(ctypes.c_double),
ndpointer(ctypes.c_int32),
ndpointer(ctypes.c_int32),
ndpointer(ctypes.c_int32),
ndpointer(ctypes.c_int32),
ndpointer(ctypes.c_int32),
ndpointer(ctypes.c_int32)];


#void getAij_matrix1(double*AIJ, double*poi, int n, double*pts,int npt, int ptdim_max,int*ptdim, double*ORIGS, double*csysesLoc,double*D,double*S,int*gid,int*isFirst,int*isLast,int*isWake){
#void getAij_matrix1(double*AIJ, double*poi, int m, double*pts, double*ORIGS, double*csysesLoc, int nw, double*ptsW, double*ORIGSW, double*csysesLocW,int*firstPanel,int*lastPanel){

testlib.getAij_matrix1.restype = ctypes.c_void_p;
testlib.getAij_matrix1.argtypes = [ndpointer(ctypes.c_double),
ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
ctypes.c_int32,
ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
ctypes.c_int32,
ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
ndpointer(ctypes.c_int32),
ndpointer(ctypes.c_int32),
ndpointer(ctypes.c_int32),
ndpointer(ctypes.c_int32),
ndpointer(ctypes.c_int32)];


#def getBoundingBoxes(*args):
	#bbs=[]
	#for wing in args:
		#bb=[getBoundingBox(wing)[0],getBoundingBox(wing)[1]]
		#bbs.append(bb)
	#return bbs

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


class Wake1(object):
	def __init__(self):
		self.dicts={}
		self.rowsMax=0
		self.restricted=np.array(())
		pass
	
	def addDict(self,d):
		self.dicts[len(self.dicts)]=d
		if len(d) > self.rowsMax:
			self.rowsMax=len(d)
	
	def translate(self,pr1,v):
		for i in range(0,self.rowsMax):
			for w in self.dicts:
				d=self.dicts[w]
				dist=d[i]['displDistribution']
				for j in range(0,len(d[i]['rowPanels'])):
					p=pr1.getPanelByGid(d[i]['rowPanels'][j])
					for l in range(0,4):
						if not pointInArray(self.restricted,p.pts[l]):
							p.pts[l]+=v
					#for k in dd:
						#pp=pr1.getPanelByGid(k)
						#if not pointInArray(self.restricted,pp.pts[dist[0][0]]):
							#pp.pts[dist[0][0]]+=v
						##else:
							##print pp.pts[dist[0][0]]
						#if not pointInArray(self.restricted,pp.pts[dist[0][1]]):
							#pp.pts[dist[0][1]]+=v
						##else:
							##print pp.pts[dist[0][1]]
						#if not pointInArray(self.restricted,pp.pts[dist[1][0]]):
							#pp.pts[dist[1][0]]+=v
						##else:
							##print pp.pts[dist[1][0]]
						#if not pointInArray(self.restricted,pp.pts[dist[1][1]]):
							#pp.pts[dist[1][1]]+=v
						##else:
							##print pp.pts[dist[1][1]]









		
	def findRowDisplacement(self,pr1,i):
		for w in self.dicts:
			d=self.dicts[w]
			#try:
			#d[i]
			poiRow=np.zeros((len(d[i]['rowPanels'])*2,3))
			for j in range(0,len(d[i]['rowPanels'])):
				p=pr1.getPanelByGid(d[i]['rowPanels'][j])
				poiRow[j*2+0]=p.pts[d[i]['pids'][0]]
				poiRow[j*2+1]=p.pts[d[i]['pids'][1]]
				
			vRow=pr1.getVelocityAtPoints(poiRow)+np.cross(pr1.omega,poiRow-pr1.COR)

			for j in range(0,len(d[i]['rowPanels'])):
				#d[i]['rowDisplacement'][j]={}
				p=pr1.getPanelByGid(d[i]['rowPanels'][j])
				poi=p.pts[d[i]['pids'][0]]
				vv=vRow[j*2+0]
				#vv=pr1.getVelocityAtPoint(poi)
				dx1=p.pts[2,0]-p.pts[1,0]
				dy1=p.pts[2,1]-p.pts[1,1]
				dz1=p.pts[2,2]-p.pts[1,2]
				v=np.array([0,vv[1]/vv[0]*dx1,vv[2]/vv[0]*dx1])
				
				
				d[i]['rowDisplacement'][j][0]=v-np.array([0,dy1,dz1])

				dx2=p.pts[3,0]-p.pts[0,0]
				dy2=p.pts[3,1]-p.pts[0,1]
				dz2=p.pts[3,2]-p.pts[0,2]
				
				
				poi=p.pts[d[i]['pids'][1]]
				vv=vRow[j*2+1]
				#vv=pr1.getVelocityAtPoint(poi)
				v=np.array([0,vv[1]/vv[0]*dx2,vv[2]/vv[0]*dx2])
				d[i]['rowDisplacement'][j][1]=v-np.array([0,dy2,dz2])
			#except:
				#pass
	def allignRowToFreestream(self,pr1,i,wv=None):
		for w in self.dicts:
			d=self.dicts[w]
			try:
				#d[i]
					
				for j in range(0,len(d[i]['rowPanels'])):
					p=pr1.getPanelByGid(d[i]['rowPanels'][j])
					if wv is not None:
						vv=wv
					else:
						vv=pr1.getFreeVelocity()
					dx1=p.pts[2,0]-p.pts[1,0]
					dy1=p.pts[2,1]-p.pts[1,1]
					dz1=p.pts[2,2]-p.pts[1,2]
					v=np.array([0,vv[1]/vv[0]*dx1,vv[2]/vv[0]*dx1])
					
					
					d[i]['rowDisplacement'][j][0]=v-np.array([0,dy1,dz1])
					#d[i]['rowDisplacement'][j][0]=v-np.array([0,0,dz1])
					
					dx2=p.pts[3,0]-p.pts[0,0]
					dy2=p.pts[3,1]-p.pts[0,1]
					dz2=p.pts[3,2]-p.pts[0,2]
					
					
					v=np.array([0,vv[1]/vv[0]*dx2,vv[2]/vv[0]*dx2])
					d[i]['rowDisplacement'][j][1]=v-np.array([0,dy2,dz2])
					#d[i]['rowDisplacement'][j][1]=v-np.array([0,0,dz2])

			except:
				print w,i,j#, d
				
	def setHorseshoeWake(self,pr1):
		for i in self.dicts:
			for j in self.dicts[i]:
				for k in self.dicts[i][j]['rowPanels']:
					pid=self.dicts[i][j]['rowPanels'][k]
					p=pr1.getPanelByGid(pid)
					p.activeEdges[0]=0
					p.activeEdges[2]=0
			for k in self.dicts[i][0]['rowPanels']:
				pid=self.dicts[i][0]['rowPanels'][k]
				p=pr1.getPanelByGid(pid)
				p.activeEdges[0]=1


	def deformRow(self,pr1,i,kw,freeDOF=np.array([1,1,1])):
		for w in self.dicts:
			d=self.dicts[w]
			#try:
			d[i]['rowDisplacement']
			dist=d[i]['displDistribution']
			for j in range(0,len(d[i]['rowPanels'])):
				p=pr1.getPanelByGid(d[i]['rowPanels'][j])
				if not pointInArray(self.restricted,p.pts[d[i]['pids'][0]]):
					p.pts[d[i]['pids'][0]]+=np.multiply(d[i]['rowDisplacement'][j][0]*kw,freeDOF)
				if not pointInArray(self.restricted,p.pts[d[i]['pids'][1]]):	
					p.pts[d[i]['pids'][1]]+=np.multiply(d[i]['rowDisplacement'][j][1]*kw,freeDOF)
				dd=d[i]['nextCol'][j]
				
				for k in dd:
					pp=pr1.getPanelByGid(k)
					if not pointInArray(self.restricted,pp.pts[dist[0][0]]):
						pp.pts[dist[0][0]]+=np.multiply(d[i]['rowDisplacement'][j][0]*kw,freeDOF)
					#else:
						#print pp.pts[dist[0][0]]
					if not pointInArray(self.restricted,pp.pts[dist[0][1]]):
						pp.pts[dist[0][1]]+=np.multiply(d[i]['rowDisplacement'][j][0]*kw,freeDOF)
					#else:
						#print pp.pts[dist[0][1]]
					if not pointInArray(self.restricted,pp.pts[dist[1][0]]):
						pp.pts[dist[1][0]]+=np.multiply(d[i]['rowDisplacement'][j][1]*kw,freeDOF)
					#else:
						#print pp.pts[dist[1][0]]
					if not pointInArray(self.restricted,pp.pts[dist[1][1]]):
						pp.pts[dist[1][1]]+=np.multiply(d[i]['rowDisplacement'][j][1]*kw,freeDOF)
					#else:
						#print pp.pts[dist[1][1]]


			#except:
				#pass


	def relax(self,pr1,kw=0.1,findDiplacement=True,deform=True,freeDOF=np.array([1,1,1])):
		for i in range(0,self.rowsMax):
			if findDiplacement:
				self.findRowDisplacement(pr1,i)
			if deform:
				self.deformRow(pr1,i,kw,freeDOF)
			
	def allignToFreestream(self,pr1,wv=None,freeDOF=np.array([1,1,1])):
		for i in range(0,self.rowsMax):
			self.allignRowToFreestream(pr1,i,wv)
			self.deformRow(pr1,i,1.,freeDOF)



class subSection(GlobalEnvironment):
	def __init__(self):
		self.nSpanwise=0
		self.nPanChordwise=0
		self.nLPanChordwise=0
		self.nPanels=0
		self.panels={}
		self.nLPanels=0
		self.LPanels={}
		self.refPT=np.zeros(3)
		self.blIter=0

	def addPanel(self,panel):
		self.panels[self.nPanels]=panel
		self.nPanels=self.nPanels+1

		if not panel.isWake:
			self.LPanels[self.nLPanels]=panel
			self.nLPanels+=1;

	def computePlainGridParams(self):
		n=self.nLPanels
		p0=self.panels[0]
		for i in range(2,n):
			if commonEdgeTwoPanels(p0,self.panels[i]):
				break
		self.nLPanChordwise=i+1
		if self.nLPanels%self.nLPanChordwise !=0:
			p0=self.panels[n-1]
			j=3
			for i in range(n-3,-1,-1):
				if commonEdgeTwoPanels(p0,self.panels[i]):
					break
				else:
					j+=1
			self.nLPanChordwise=j#+1			
		self.nSpanwise=self.nLPanels/self.nLPanChordwise
	
	def joinChordwise(self,other,reverseSpanwise=False,reverseChordwiseSelf=False,reverseChordwiseOther=False):
		LPanels={}
		panels={}
		nSpanwise=self.nSpanwise
		rngSpan=range(0,nSpanwise)
		if reverseSpanwise:
			rngSpan=range(nSpanwise-1,-1,-1)
		rngChord=range(0,self.nLPanChordwise)
		rngChordOther=range(0,other.nLPanChordwise)
		if reverseChordwiseOther:
			rngChordOther=range(other.nLPanChordwise-1,-1,-1)	
		if reverseChordwiseSelf:
			rngChord=range(self.nLPanChordwise-1,-1,-1)
			
		pc=0
		for i in rngSpan:
			for j in  rngChord:
				LPanels[pc]=copy.copy(self.LPanels[self.nLPanChordwise*i +j])
				pc+=1
			for j in  rngChordOther:
				LPanels[pc]=copy.copy(other.LPanels[other.nLPanChordwise*i +j])
				pc+=1
		self.__init__()
		for i in LPanels:
			self.addPanel(copy.copy(LPanels[i]))
		self.nLPanChordwise=self.nLPanChordwise+other.nLPanChordwise
		self.nSpanwise=nSpanwise
		
	def setPGIDS(self):
		for i in range(0,len(self.panels)):
			self.panels[i].pgid=i
	
	def getPanelByPGID(self,i):
		return self.panels[i]

	def getGidsInsidePointCloud(self,points,addWakes=False):
		def pnt_in_cvex_hull_1(hull, pnt):
			'''
			Checks if `pnt` is inside the convex hull.
			`hull` -- a QHull ConvexHull object
			`pnt` -- point array of shape (3,)
			'''
			new_hull = ConvexHull(np.concatenate((hull.points, [pnt])))
			if np.array_equal(new_hull.vertices, hull.vertices): 
				return True
			return False


		hull = ConvexHull(points, incremental=True)
		hull_points = hull.points[hull.vertices, :]
		gids=[]
		for i in self.panels:
			if pnt_in_cvex_hull_1(hull,self.panels[i].getCpoint()):
				if not addWakes:
					gids.append(self.panels[i].gid)
		return gids

	def getChord(self):
		p0=self.getLiftPanel(0).getCpoint()
		dv=np.zeros(self.getNumberOfLiftPanels())
		for i in range(0,self.getNumberOfLiftPanels()):
			dv[i]=np.linalg.norm(self.getLiftPanel(i).getCpoint()-p0)
		return dv.max()
	
	def sortPanels(self,dom1,mode='addedWakes'):
		strips=self.getSpanwiseStrips(liftOnly=True)
		panels={}
		ipan=0
		for i in strips:
			s=strips[i]
			for j in s.panels:
				p=s.panels[j]
				panels[ipan]=p
				ipan+=1
			if p.isLast:
				for k in p.wakes:
					panels[ipan]=dom1.getPanelByGid(k)
					ipan+=1
			else:
				print('sortPanels failed: bad panels order')
		self.panels=panels	
		self.nPanChordwise=self.nLPanChordwise+len(p.wakes)


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



	def getBLStrip(self,dom1,sdir='up',reverse=False):
		#print 'extracting bl strip'
		vref=np.linalg.norm(dom1.getFreeVelocity())
		n=self.getNumberOfPanels()
		nl=self.getNumberOfLiftPanels()

		VRef=np.linalg.norm(dom1.getFreeVelocity())
		
		pp0=[]
		cps=np.zeros(n)
		v=np.zeros(n)
		
		for i in range(0,n):
			p=self.getPanel(i)
			cps[i]=p.getPressureCoef()

			pp0.append(p)
			if not p.isWake:
				if dom1.blIter>0:# and p.displThick<0.1:
					p.uebl=np.copy(np.linalg.norm(p.vloc[[0,1]]))/VRef
				else:
					p.uebl=np.copy(np.linalg.norm(p.vloc[[0,1]]))/VRef

			if p.isWake:
				if sdir=='up':
					#p.ue=np.linalg.norm(dom1.getVelocityAtPoints(p.getCpoint()- p.normalVect*p.blDtUp))
					p.uebl=p.blUeUp/VRef
				else:
					#p.ue=np.linalg.norm(dom1.getVelocityAtPoints(p.getCpoint()+ p.normalVect*p.blDtDown))
					p.uebl=p.blUeDown/VRef
			#p.w=0.0
			#p.uebl=np.copy(p.uebl/VRef)
			v[i]=p.uebl
		n04=nl/4
		iStag=v[n04:nl-n04].argmin()+n04
		#print iStag
		#self.v=v
		#self.iStag=iStag
		idxx=np.zeros(n,dtype=np.bool)
		
		if sdir=='down':
			idxx[iStag-1::]=True
		else:
			idxx[0:iStag]=True
			idxx[nl::]=True
		panels=[]
		for i in range(0,n):
			if idxx[i]:
				panels.append(pp0[i])
				
		npan=len(panels)
		nls=0
		for i in range(0,npan):
			p=panels[i]
			if not p.isWake:
				nls+=1	
				
		if sdir=='up':
			pps0=pp0[iStag+1]
		else:
			pps0=pp0[iStag-2]
			
		if reverse:
			panels[0:nls]=panels[0:nls][::-1]
		pts=np.zeros((npan+1,2))
		uebl=np.zeros(npan+1)
		dlspbl=np.zeros(npan+1)
		vnpbl=np.zeros(npan+1)
		uedispbl=np.zeros(npan+1)
		cfsbl=np.zeros(npan+1)
		
		for i in range(0,npan):
			p=panels[i]
			#p.w=0.0
			#p.displThick=0.0

			pts[i+1]=np.copy(p.getCpoint()[[0,2]])
			uebl[i+1]=np.copy(p.uebl)
			dlspbl[i+1]=np.copy(p.dlspbl)
			vnpbl[i+1]=np.copy(p.vnpbl)
			uedispbl[i+1]=np.copy(p.ue)#dispbl)
			cfsbl[i+1]=np.copy(p.cfsbl)

		
		pts[0]=np.copy((pps0.getCpoint()[[0,2]]-panels[0].getCpoint()[[0,2]])*0.5+ panels[0].getCpoint()[[0,2]])
		uebl[0]=(pps0.uebl+panels[0].uebl)*0.5
		dlspbl[0]=(pps0.dlspbl+panels[0].dlspbl)*0.5
		vnpbl[0]=(pps0.vnpbl+panels[0].vnpbl)*0.5
		uedispbl[0]=(pps0.uedispbl+panels[0].uedispbl)*0.5
		cfsbl[0]=pps0.cfsbl
		
		s=IblmStrip()
		if sdir=='up':
			s.sdir='up'
		else:
			s.sdir='down'
		pts[:,0]-=pts[0,0]#pts[:,0].min()
		pts[:,1]-=pts[0,1] #pts[:,1].min()
		pts/=self.getChord()
		s.VRef=np.linalg.norm(dom1.getFreeVelocity())
		s.LRef=self.getChord()
		s.panels=panels
		s.setPts(pts)
		s.setUebl(uebl)
		s.setNlpt(nls)
		s.setArray('dlspbl',dlspbl)
		s.setArray('vnpbl',vnpbl)
		s.setArray('uedispbl',uedispbl)
		s.setArray('cfsbl',cfsbl)
		
		s.irestart= self.blIter


		
		return s


	
	def examineBLAnalitic(self,dom1):
		vref=np.linalg.norm(dom1.getFreeVelocity())

		blStripsUp={}
		blStripsDown={}
		self.blStripsDown=blStripsDown
		self.blStripsUp=blStripsUp

		secs=self.getSpanwiseStrips()
		self.secs=secs
		rho=secs[0].getLiftPanel(0).getRho()
		mu=secs[0].getLiftPanel(0).getMu()
		
		for i in secs:
			s=secs[i]
			
			blStripsUp[i]=s.getBLStrip(dom1,'up',True)
			blStripsUp[i].VRef=vref
			blStripsUp[i].rho=rho
			blStripsUp[i].mu=mu
			blStripsUp[i].setBLParam(dom1)
			
			blStripsDown[i]=s.getBLStrip(dom1,'down',False)
			blStripsDown[i].VRef=vref
			blStripsDown[i].rho=rho
			blStripsDown[i].mu=mu
			blStripsDown[i].setBLParam(dom1)
		
	def setBlowing(self):
		#for i in self.blStripsUp:
			#self.blStripsUp[i].setBlowing()
		#for i in self.blStripsDown:
			#self.blStripsDown[i].setBlowing()
		pass
	
	def setRefPT(self,PT):
		self.refPT=PT

	def getPanelList(self):
		return self.panels.keys()

	def getNumberOfPanels(self):
		return len(self.panels)

	def getPanel(self,i):
		return self.panels[i]

	def getLiftPanelList(self):
		return self.LPanels.keys()

	def getNumberOfLiftPanels(self):
		return len(self.LPanels)

	def getLiftPanel(self,i):
		return self.LPanels[i]
	
	def translate(self,v):
		for i in self.getPanelList():
			self.getPanel(i).pts+=v
			
	def getBoundingBox(self,scl=1.0,liftOnly=True):
		if liftOnly:
			points=np.zeros((len(self.getLiftPanelList())*4,3))
			for i in self.getLiftPanelList():
				points[i*4:(i+1)*4]=self.getLiftPanel(i).pts
		else:
			points=np.zeros((len(self.getPanelList())*4,3))
			for i in self.getPanelList():
				points[i*4]=self.getPanel(i).pts
		#a = zeros((2,2))
		#a[:,0] = np.min(points, axis=0)
		#a[:,1] = np.max(points, axis=0)
		A = np.min(points, axis=0)
		B = np.max(points, axis=0)
		v=A-B
		return A-v*scl,B+v*scl

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
		for i in self.getLiftPanelList():
			if not self.getPanel(i).isWake:
				A+=self.getPanel(i).getArea()
		return A		

	def getPanelGidsByBoundingBox(self,A,B,addWakes=True,onlyWakes=True):
		l=[]
		for i in self.getPanelList():
			panelI=self.getPanel(i)
			p=panelI#.getCpoint()
			panelI=self.getPanel(i)
			c1=True
			c2=True
			for j in range(0,3):
				c1*=np.any(p.pts[:,j]>=A[j]) and np.any(p.pts[:,j]<=B[j])
				c2*=np.any(p.pts[:,j]<=A[j]) and np.any(p.pts[:,j]>=B[j])
			#if (p[0] > A[0] and p[0] < B[0]) and (p[1] > A[1] and p[1] < B[1]) and (p[2] > A[2] and p[2] < B[2]):
			if c1 or c2:
				if not panelI.isWake and not onlyWakes:
					l.append(panelI.gid)
				if addWakes:
					l.append(panelI.gid)
		return l

	def getPanelsByBoundingBox(self,A,B,addWakes=True,onlyWakes=True):
		l=PanelGroup()
		for i in self.getPanelList():
			panelI=self.getPanel(i)
			pt=panelI.getCpoint()
			c1=True
			c2=True
			for j in range(0,3):
				c1*=pt[j]>A[j] and pt[j]<B[j]
				c2*=pt[j]<A[j] and pt[j]>B[j]
			#if (p[0] > A[0] and p[0] < B[0]) and (p[1] > A[1] and p[1] < B[1]) and (p[2] > A[2] and p[2] < B[2]):
			if c1 or c2:
				if not panelI.isWake and not onlyWakes:
					l.addPanel(panelI)
				if addWakes:
					l.addPanel(panelI)
		return l


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

	def getWakeDict(self,globalRowOffset=0):
		nr=self.nPanChordwise - self.nLPanChordwise
		d={}
		for i in range(0,nr):
			d[i]={}
			d[i]['rowPanels']={}
			d[i]['nextCol']={}
			d[i]['pids']=(2,3)
			d[i]['displDistribution']=( (1,2),(0,3))
			d[i]['rowDisplacement']={}
			d[i]['vOld']={}
			for j in range(0,self.nSpanwise):
				d[i]['rowPanels'][j]=self.panels[self.nPanChordwise*j +self.nLPanChordwise + i].gid
				d[i]['nextCol'][j]=[]
				for k in range(i+1,nr):
					d[i]['nextCol'][j].append(self.panels[self.nPanChordwise*j +self.nLPanChordwise + k].gid)
				
				d[i]['rowDisplacement'][j]={}
				d[i]['rowDisplacement'][j][0]=np.zeros(3)
				d[i]['rowDisplacement'][j][1]=np.zeros(3)
		
				d[i]['vOld'][j]={}
				d[i]['vOld'][j][0]=np.zeros(3)
				d[i]['vOld'][j][1]=np.zeros(3)

		return d
	
	def getPTS2Array(self):
		r=np.zeros((self.nLPanels*4,3))
		for i in range(0,self.nLPanels):
			p=self.getLiftPanel(i)
			for j in range(0,4):
				r[i*4+j]=p.pts[j]
		return r

	def getBodyPanels(self):
		return self.LPanels

class PanelGroup(subSection):
	def addSubsection(self,r,join=True):
		for i in r.getPanelList():
			self.addPanel(r.getPanel(i))
		if join:
			self.nSpanwise+=r.nSpanwise
			self.nLPanChordwise=r.nLPanChordwise
		
	def addPanelGroup(self,pg):
		self.addSubsection(pg)

class GlobalDomain(GlobalEnvironment):
	def __init__(self):
		self.nRuns=0
		self.view=PanelViewer3D()
		self.refPT=np.zeros(3)
		self.regions={}
		self.gids={}
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
				if self.getRegion(r).getPanel(i).isWake:
					self.wakePanels.append(gid) 
		
	def getPanelByGid(self,gid):
		return self.gids[gid]
				
	def setRefPT(self,PT):
		self.refPT=PT

	def addRegion(self,reg):
		self.regions[self.nReg]=reg
		self.nPan=self.nPan+self.regions[self.nReg].getNumberOfPanels()
		self.nLPan+=self.regions[self.nReg].getNumberOfLiftPanels()
		
		self.buildGids()
		self.nReg=self.nReg+1
#		print '===================== add region ================================='
#		print reg 
#		print self.regions[self.nReg]
#
#		for i in range(0,len(reg.panels)):
#			print reg.panels[i]
#			print self.regions[self.nReg].panels[i]
		
	def getForceAtPTByGidList(self,PT,l):
		f=np.zeros(3)
		m=np.zeros(3)
		for i in self.gids.keys():
			panelI=self.getPanelByGid(i)
			if panelI.gid in l:
				f=f+panelI.getForce()
				r=panelI.midPoint -PT
				m=m+np.cross(r,panelI.getForce())
		return np.hstack((f,m))

	def getPanelGidsByBoundingBox(self,A,B,addWakes=True,onlyWakes=True):
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
					if not p.isWake and not onlyWakes:
						l.append(p.gid)
					if addWakes:
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



	def getPHI0AtPoint(self,PT):
		return np.dot(self.getFreeVelocity(),PT)



	def getPHIAtPoint(self,PT):
		phi=self.getPHI0AtPoint(PT)
	#	phi=0.0
		for r in self.getRegionList():
			for i in self.getRegion(r).getPanelList():
				phi=phi+self.getRegion(r).getPanel(i).getPHI(PT)
		return phi
	
	def getPHIInternalAtPoint(self,PT):
		phi=0.0
		for r in self.getRegionList():
			for i in self.getRegion(r).getPanelList():
				phi=phi+self.getRegion(r).getPanel(i).getPHI(PT)
		return phi

		

	def reinitPanelsGeometry(self):
		for r in self.getRegionList():
			for i in self.getRegion(r).getPanelList():
				self.getRegion(r).getPanel(i).initGeometry()

	def setPanelsG(self,G):
		for r in self.getRegionList():
			for i in self.getRegion(r).getPanelList():
				self.getRegion(r).getPanel(i).setG(G)

	def setPanelsFreeStream(self,v):
		for r in self.getRegionList():
			for i in self.getRegion(r).getPanelList():
				self.getRegion(r).getPanel(i).setFreeStream(v)

	def getPanelOffset(self,i):
		ofs=0
		for i in range(0,i):
			ofs=ofs+self.getRegion(i).getNumberOfPanels()
		return ofs

	def getLiftPanelOffset(self,i):
		ofs=0
		for i in range(0,i):
			ofs=ofs+self.getRegion(i).getNumberOfLiftPanels()
		return ofs


	def setWakeRC(self):#,RC):
		#self.RCAll=np.ones(self.nPan,dtype=np.float64)
		v=np.linalg.norm(self.getFreeVelocity())
		gids=self.gids.keys()
		for i in range(0,self.nPan):	
			p=self.getPanelByGid(gids[i])
			if p.isWake:
				#pf=self.getPanelByGid(p.firstPanel)
				#x=np.linalg.norm(p.getCpoint()-pf.getCpoint())
				#t=x/v
				#t=1.
				#self.RCAll[i]=np.sqrt(4.*self.nu*t)
				self.RCAll[i]=self.vortexRC

	def allocPointsAll(self):

		self.ptsAll=np.zeros((self.nPan,3*4),dtype=np.float64)

		self.ptdimAll=np.ones(self.nPan,dtype=np.int32)*4
		self.origsAll=np.zeros((self.nPan,3),dtype=np.float64)
		self.csysesLocAll=np.zeros((self.nPan,9),dtype=np.float64)

		self.isWakeAll= np.ones(self.nPan,dtype=np.int32)*-1

		self.DAll=np.zeros(self.nPan,dtype=np.float64)
		self.SAll=np.zeros(self.nPan,dtype=np.float64)
		#self.RCAll=np.zeros(self.nPan,dtype=np.float64)
		self.RCAll=np.ones(self.nPan,dtype=np.float64)
		self.setWakeRC()

		self.activeEdgeAll=np.ones((self.nPan,4),dtype=np.int32)
		gids=self.gids.keys()
		for i in range(0,self.nPan):	
			p=self.getPanelByGid(gids[i])
			self.ptsAll[i]=np.float64( p.pts.flatten())
			self.origsAll[i]=np.float64(p.controlPoint)
			self.csysesLocAll[i]=np.float64(p.csys.flatten())
			self.SAll[i]=np.float64(p.S)
			self.DAll[i]=np.float64(p.G)
			self.activeEdgeAll[i]=p.activeEdges

	def solveAij00(self,normalsOrientation=1.0):
		self.buildGids()
		print 'assembling Aij matrix'
		print strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())

####
####		for r in self.getRegionList():
####			nPan=self.getRegion(r).getNumberOfLiftPanels()
####			ofs=self.getLiftPanelOffset(r)	
####			for i in range(0,nPan):
####				p=self.getRegion(r).getLiftPanel(i)
#####				self.getRegion(r).getLiftPanel(i).setS(normalsOrientation*np.dot(p.normalVect,self.getFreeVelocity()))
####				vloc=p.transform2loc(np.zeros(3),self.getFreeVelocity())
####				self.getRegion(r).getLiftPanel(i).setS(normalsOrientation * vloc[2])
####		for r1 in self.getRegionList():
####			for r2 in self.getRegionList():
####				ofsI=self.getPanelOffset(r1)
####				nPanI=self.getRegion(r1).getNumberOfPanels()
####				for i in range(0,nPanI):
####					panelI=self.getRegion(r1).getPanel(i)
####					panelI.G=0.0
####
####
####		A=np.zeros((self.nPan,self.nPan))
####		B=np.zeros(self.nPan)
####
####		for r1 in self.getRegionList():
####			for r2 in self.getRegionList():
####				ofsI=self.getPanelOffset(r1)
####				nPanI=self.getRegion(r1).getNumberOfPanels()
####				for i in range(0,nPanI):
####					panelI=self.getRegion(r1).getPanel(i)
####					B[i+ofsI]=self.getPHIInternalAtPoint(panelI.controlPoint)
####
####
####		for r1 in self.getRegionList():
####			for r2 in self.getRegionList():
####				ofsI=self.getPanelOffset(r1)
####				ofsJ=self.getPanelOffset(r2)
####				nPanI=self.getRegion(r1).getNumberOfPanels()
####				nPanJ=self.getRegion(r2).getNumberOfPanels()
####				for i in range(0,nPanI):
####					panelI=self.getRegion(r1).getPanel(i)
####					for j in range(0,nPanJ):
####						panelJ=self.getRegion(r2).getPanel(j)
####						if not panelI.isWake:
####							A[i+ofsI][j+ofsJ]=panelJ.getDoubletPHI(panelI.getCpoint(),1.0) 
####						else:
####							if panelJ.isWake:
####								A[j+ofsJ][j+ofsJ]=1.0
####								
####							if panelJ.isFirst:
####								if panelI.firstPanel == panelJ.gid:
####									A[i+ofsI][j+ofsJ]=1.0
####
####							if panelJ.isLast:
####								if panelI.lastPanel == panelJ.gid:
####									A[i+ofsI][j+ofsJ]=-1.0
####
####					if panelI.isWake:
####						B[i+ofsI]=0.0
####
#######################################################################################################################################				
####		self.B=-B#(1.)*np.dot(self.b,self.d)
#####		self.B= -np.dot(self.s,self.b)
####		self.B0=B
####		self.A=A
####		print 'Aij matrix assembled'
####		print 'solving Aij'
####
####		self.G=lgmres(self.A,self.B,tol=1e-16)[0]
#####		self.G=np.squeeze(np.linalg.solve(A,self.B))
#####		self.G=np.squeeze(np.linalg.lstsq(A,self.B)[0])
####		print 'postprocessing results'
####		print strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
####
#####		for r in self.getRegionList():
#####			nPan=self.getRegion(r).getNumberOfLiftPanels()
#####			ofs=self.getLiftPanelOffset(r)	
#####			for i in range(0,nPan):
#####				self.getRegion(r).getPanel(i).setG(self.G[ofs+i])
#####				p=self.getRegion(r).getLiftPanel(i)
#####				if p.isLast:
#####					if p.wakes:
#####						g= p.G - self.getPanelByGid(p.firstPanel).G #- p.G
#####						for ppw in range(0,len(p.wakes)):
#####							self.getPanelByGid(p.wakes[ppw]).setG(g)
#####
####		for r in self.getRegionList():
####			nPan=self.getRegion(r).getNumberOfPanels()
####			ofs=self.getPanelOffset(r)	
####			for i in range(0,nPan):
####				self.getRegion(r).getPanel(i).setG(self.G[ofs+i])
####				#########################
####
####
####
####		print 'panels G set'
####		print strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())			
####		for r in self.getRegionList():
####			nPan=self.getRegion(r).getNumberOfPanels()
####			ofs=self.getPanelOffset(r)
####			for i in range(0,nPan):	
####				panelI=self.getRegion(r).getPanel(i)
####				CPT=panelI.getCpoint()
####		#		panelI.phiInt=self.getPHIInternalAtPoint(CPT)
#####				panelI.setFreeStream(self.getVelocityAtPoint(CPT))
####		for r in self.getRegionList():
####			nPan=self.getRegion(r).getNumberOfPanels()
####			ofs=self.getPanelOffset(r)
####			for i in range(0,nPan):	
####				panelI=self.getRegion(r).getPanel(i)
####				panelI.setForce()
####		print 'results computed'
####		print strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
###### 
#void getPHI_at_points(double*PHI, double*poi, int n, double*pts,int npt, int ptdim_max,int*ptdim, double*ORIGS, double*csysesLoc,double*D,double*S){
#		testlib.getPHI_at_points.restype = ctypes.c_void_p;
#		testlib.getPHI_at_points.argtypes = [ndpointer(ctypes.c_double),
#		ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
#		ctypes.c_int32,
#		ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
#		ctypes.c_int32,
#		ctypes.c_int32,
#		ndpointer(ctypes.c_int32),
#		ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
#		ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
#		ndpointer(ctypes.c_double),
#		ndpointer(ctypes.c_double)];	

		self.B=np.zeros(self.nPan,dtype=np.float64)
		poi=np.zeros((self.nPan,3),dtype=np.float64)
		n=self.nPan
		pts=np.zeros((self.nPan,3*4),dtype=np.float64)
		npt=self.nPan
		ptdim_max=4
		ptdim=np.ones(self.nPan,dtype=np.int32)*4
		origs=np.zeros((self.nPan,3),dtype=np.float64)
		csysesLoc=np.zeros((self.nPan,9),dtype=np.float64)
		D=np.zeros(self.nPan)
		SO=np.zeros(self.nPan)
		print 'assembling RHS by testlib '+strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
		for r in self.getRegionList():
			nPan=self.getRegion(r).getNumberOfPanels()
			ofsI=self.getPanelOffset(r)
			for i in range(0,nPan):	
				p=self.getRegion(r).getPanel(i)
				poi[i+ofsI]=np.float64(p.controlPoint)
				pts[i+ofsI]=np.float64( p.pts.flatten())
				origs[i+ofsI]=np.float64(p.controlPoint)
				csysesLoc[i+ofsI]=np.float64(p.csys.flatten())
				D[i+ofsI]=np.float64(0.0)
				#vloc=p.transform2loc(np.zeros(3),self.getFreeVelocity())
				if p.isWake:
					SO[i+ofsI]=np.float64(0.0-p.getBlowingVelocity())
				else:
					rvector=p.getCpoint()-self.COR
					vrot=np.cross(self.omega,rvector)
					#SO[i+ofsI]=np.float64((normalsOrientation * vloc[2])  -p.getBlowingVelocity() )
					SO[i+ofsI]=np.float64(np.dot((self.getFreeVelocity()+vrot) ,p.normalVect))
		with nogil:
			testlib.getPHI_at_points(self.B,poi,n,pts,npt,ptdim_max,ptdim,origs,csysesLoc,D,SO)
		for r in self.getRegionList():
			nPan=self.getRegion(r).getNumberOfPanels()
			ofsI=self.getPanelOffset(r)
			for i in range(0,nPan):	
				p=self.getRegion(r).getPanel(i)
				if p.isWake:
					self.B[i+ofsI]=0.0	
		self.SO=SO

		self.A=np.zeros(self.nPan*self.nPan ,dtype=np.float64)
		poi=np.zeros((self.nPan,3),dtype=np.float64)
		n=self.nPan#*self.nPan
		pts=np.zeros((self.nPan,3*4),dtype=np.float64)
		npt=self.nPan
		ptdim_max=4
		ptdim=np.ones(self.nPan,dtype=np.int32)*4
		origs=np.zeros((self.nPan,3),dtype=np.float64)
		csysesLoc=np.zeros((self.nPan,9),dtype=np.float64)
		D=np.zeros(self.nPan,dtype=np.float64)
		S=np.zeros(self.nPan,dtype=np.float64)
		gid=np.zeros(self.nPan,dtype=np.int32)
		isFirst=np.ones(self.nPan,dtype=np.int32)*-1
		isLast= np.ones(self.nPan,dtype=np.int32)*-1
		isWake= np.ones(self.nPan,dtype=np.int32)*-1
		firstPanel=np.ones(self.nPan,dtype=np.int32)*-1
		lastPanel=np.ones(self.nPan,dtype=np.int32)*-1
		print 'assembling AIJ by testlib '+strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
		for r in self.getRegionList():
			nPan=self.getRegion(r).getNumberOfPanels()
			ofsI=self.getPanelOffset(r)
			for i in range(0,nPan):	
				p=self.getRegion(r).getPanel(i)

				poi[i+ofsI]=np.float64(p.controlPoint)
				pts[i+ofsI]=np.float64( p.pts.flatten())
				origs[i+ofsI]=np.float64(p.controlPoint)
				csysesLoc[i+ofsI]=np.float64(p.csys.flatten())
				D[i+ofsI]=np.float64(1.0)
				S[i+ofsI]=np.float64(0.0)

				gid[i+ofsI]=np.int32(p.gid)
				if p.isFirst:
					isFirst[i+ofsI]=np.int32(1)
				if p.isLast:
					isLast[i+ofsI]=np.int32(1)
				if p.isWake:
					isWake[i+ofsI]=np.int32(1)
				if p.firstPanel != None:
					firstPanel[i+ofsI]=np.int32(p.firstPanel)
				if p.lastPanel != None:
					lastPanel[i+ofsI]=np.int32(p.lastPanel)

		testlib.getAij_matrix(self.A,poi,n,pts,npt,ptdim_max,ptdim,origs,csysesLoc,D,S,gid,isFirst,isLast,isWake,firstPanel,lastPanel)
		del poi,pts,ptdim,origs,csysesLoc, D, S, gid,isFirst, isLast, isWake, firstPanel, lastPanel
		
		self.A=self.A.reshape((self.nPan,self.nPan))
#		self.islast=isLast
#		self.lastPanel=lastPanel
#		for r1 in self.getRegionList():
#			for r2 in self.getRegionList():
#				ofsI=self.getPanelOffset(r1)
#				ofsJ=self.getPanelOffset(r2)
#				nPanI=self.getRegion(r1).getNumberOfPanels()
#				nPanJ=self.getRegion(r2).getNumberOfPanels()
#				for i in range(0,nPanI):
#					panelI=self.getRegion(r1).getPanel(i)
#					for j in range(0,nPanJ):
#						panelJ=self.getRegion(r2).getPanel(j)
#						if  panelI.isWake:
#							if panelJ.isWake:
#								self.A[j+ofsJ][j+ofsJ]=1.0
#								
#							if panelJ.isFirst:
#								if panelI.firstPanel == panelJ.gid:
#									self.A[i+ofsI][j+ofsJ]=1.0

#							if panelJ.isLast:
#								if panelI.lastPanel == panelJ.gid:
#									self.A[i+ofsI][j+ofsJ]=-1.0

		self.B=-self.B#(1.)*np.dot(self.b,self.d)
		print 'Aij matrix assembled'
		print 'solving Aij'

		self.G=lgmres(self.A,self.B,tol=1e-9)[0]
		print 'postprocessing results'
		print strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())

		for r in self.getRegionList():
			nPan=self.getRegion(r).getNumberOfPanels()
			ofs=self.getPanelOffset(r)
			region=self.getRegion(r)
			for i in range(0,nPan):
				region.getPanel(i).setG(self.G[ofs+i])
				region.getPanel(i).setS(self.SO[ofs+i])

		for r in self.getRegionList():
			nPan=self.getRegion(r).getNumberOfPanels()
			ofs=self.getPanelOffset(r)
			region=self.getRegion(r)
			for i in range(0,nPan):	
				region.getPanel(i).setForce()

		print 'results computed: '+strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())

	def solveAij0(self,normalsOrientation=1.0):
		self.buildGids()
		print 'assembling Aij matrix'
		print strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
		#if True:#self.nRuns < 1:
		if self.nRuns < 1:
			self.B=np.zeros(self.nPan,dtype=np.float64)
			self.poi=np.zeros((self.nPan,3),dtype=np.float64)

			self.pts=np.zeros((self.nPan,3*4),dtype=np.float64)

			self.ptdim=np.ones(self.nPan,dtype=np.int32)*4
			self.origs=np.zeros((self.nPan,3),dtype=np.float64)
			self.csysesLoc=np.zeros((self.nPan,9),dtype=np.float64)
			self.D=np.ones(self.nPan,dtype=np.float64)
			self.D0=np.zeros(self.nPan,dtype=np.float64)
			self.S=np.ones(self.nPan,dtype=np.float64)
			self.S0=np.zeros(self.nPan,dtype=np.float64)
			self.A0=np.zeros(self.nPan*self.nPan ,dtype=np.float64)
			self.A=np.zeros((self.nPan,self.nPan) ,dtype=np.float64)

			self.gid=np.zeros(self.nPan,dtype=np.int32)
			self.isFirst=np.ones(self.nPan,dtype=np.int32)*-1
			self.isLast= np.ones(self.nPan,dtype=np.int32)*-1
			self.isWake= np.ones(self.nPan,dtype=np.int32)*-1
			self.firstPanel=np.ones(self.nPan,dtype=np.int32)*-1
			self.lastPanel=np.ones(self.nPan,dtype=np.int32)*-1


		else:
			self.B.fill(0.0)
			self.poi.fill(0.0)

			self.pts.fill(0.0)

			self.ptdim.fill(4)
			self.origs.fill(0.0)
			self.csysesLoc.fill(0.0)
			self.D.fill(1.0)
			self.D0.fill(0.0)
			self.S.fill(1.0)
			self.S0.fill(0.0)
			self.A.fill(0.0)
			self.A0.fill(0.0)
			self.gid.fill(0)
			self.isFirst.fill(-1)
			self.isLast.fill(-1)
			self.isWake.fill(-1)
			self.firstPanel.fill(-1)
			self.lastPanel.fill(-1)

		self.nRuns+=1
		n=self.nPan
		npt=self.nPan
		ptdim_max=4

		self.A0.shape = (self.nPan*self.nPan)
		print 'assembling RHS by testlib '+strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
		for r in self.getRegionList():
			nPan=self.getRegion(r).getNumberOfPanels()
			ofsI=self.getPanelOffset(r)
			for i in range(0,nPan):	
				p=self.getRegion(r).getPanel(i)
				self.poi[i+ofsI]=np.float64(p.controlPoint)
				self.pts[i+ofsI]=np.float64( p.pts.flatten())
				self.origs[i+ofsI]=np.float64(p.controlPoint)
				self.csysesLoc[i+ofsI]=np.float64(p.csys.flatten())
				self.gid[i+ofsI]=np.int32(p.gid)
				#vloc=p.transform2loc(np.zeros(3),self.getFreeVelocity())
				if p.isWake:
					self.S[i+ofsI]=np.float64(0.0-p.getBlowingVelocity())
				else:
					rvector=p.getCpoint()-self.COR
					vrot=np.cross(self.omega,rvector)
					#SO[i+ofsI]=np.float64((normalsOrientation * vloc[2])  -p.getBlowingVelocity() )
					self.S[i+ofsI]=np.float64(np.dot((self.getFreeVelocity()+vrot) ,p.normalVect))

				if p.isFirst:
					self.isFirst[i+ofsI]=np.int32(1)
				if p.isLast:
					self.isLast[i+ofsI]=np.int32(1)
				if p.isWake:
					self.isWake[i+ofsI]=np.int32(1)
				if p.firstPanel != None:
					self.firstPanel[i+ofsI]=np.int32(p.firstPanel)
				if p.lastPanel != None:
					self.lastPanel[i+ofsI]=np.int32(p.lastPanel)

		affinity.set_process_affinity_mask(0,2**multiprocessing.cpu_count()-1)
		testlib.getPHI_at_points(self.B,self.poi,n,self.pts,npt,ptdim_max,self.ptdim,self.origs,self.csysesLoc,self.D0,self.S)

		print 'assembling AIJ by testlib '+strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
		affinity.set_process_affinity_mask(0,2**multiprocessing.cpu_count()-1)
		testlib.getAij_matrix0(self.A0,self.poi,n,self.pts,npt,ptdim_max,self.ptdim,self.origs,self.csysesLoc,self.D,self.S0,self.gid,self.isFirst,self.isLast,self.isWake,self.firstPanel,self.lastPanel)
		
		
		#for i in range(0,self.nPan):
			#for j in range(0,self.nPan):
				#self.A[i][j]=self.A0[i*self.nPan+j]

		self.A0.shape = (self.nPan,self.nPan)
		for r in self.getRegionList():
			nPan=self.getRegion(r).getNumberOfPanels()
			ofsI=self.getPanelOffset(r)
			for i in range(0,nPan):	
				p=self.getRegion(r).getPanel(i)
				if p.isWake:
					self.B[i+ofsI]=0.0	

		self.B*=-1.0#self.B#(1.)*np.dot(self.b,self.d)
		print 'Aij matrix assembled'
		print 'solving Aij'
		affinity.set_process_affinity_mask(0,2**multiprocessing.cpu_count()-1)
		self.G=lgmres(self.A0,self.B,tol=1e-9)[0]
		print 'postprocessing results'
		print strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())

		for r in self.getRegionList():
			nPan=self.getRegion(r).getNumberOfPanels()
			ofs=self.getPanelOffset(r)
			region=self.getRegion(r)
			for i in range(0,nPan):
				region.getPanel(i).setG(self.G[ofs+i])
				region.getPanel(i).setS(self.S[ofs+i])

		for r in self.getRegionList():
			nPan=self.getRegion(r).getNumberOfPanels()
			ofs=self.getPanelOffset(r)
			region=self.getRegion(r)
			for i in range(0,nPan):	
				region.getPanel(i).setForce()

		print 'results computed: '+strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())






	def setVelocityAtWakeBlPoints(self,wakeDistMul=100.0):
		self.allocPointsAll()
		iw=0
		print 'Solving wake bl velocity'
		for r in self.getRegionList():
			nPan=self.getRegion(r).getNumberOfPanels()
			ofsI=self.getPanelOffset(r)
			
			for i in range(0,nPan):	
				p=self.getRegion(r).getPanel(i)
				if p.isWake:
					p.blWakeSolveVelocity=False
					pFirst=self.getPanelByGid(p.firstPanel)
					if np.linalg.norm(p.getCpoint()-pFirst.getCpoint()) < np.linalg.norm(pFirst.vx)*wakeDistMul:
						p.blWakeSolveVelocity=True
						iw+=1
		#self.nWPan=iw
		self.UVW_Wake=np.zeros((iw*2,3),dtype=np.float64)
		self.poiWake=np.zeros((iw*2,3),dtype=np.float64)
		print 'Alocated '+str(iw)+' wake  points'
		npt=self.nPan
		ptdim_max=4

		
		iw=0
		for r in self.getRegionList():
			nPan=self.getRegion(r).getNumberOfPanels()
			ofsI=self.getPanelOffset(r)
			
			for i in range(0,nPan):	
				
				p=self.getRegion(r).getPanel(i)
				if p.isWake:
					if p.blWakeSolveVelocity:
						self.poiWake[iw]=np.float64(p.getCpoint()-p.getNV()*p.blDtUp)
						iw+=1
		
		#iw=self.nWPan
		for r in self.getRegionList():
			nPan=self.getRegion(r).getNumberOfPanels()
			ofsI=self.getPanelOffset(r)
			
			for i in range(0,nPan):	
				
				p=self.getRegion(r).getPanel(i)
				if p.isWake:
					if p.blWakeSolveVelocity:
						self.poiWake[iw]=np.float64(p.getCpoint()-p.getNV()*p.blDtDown)
						iw+=1
		
		print 'testlib.getUVW_at_points'
		testlib.getUVW_at_points(self.UVW_Wake,self.poiWake,iw,self.ptsAll,self.nPan,ptdim_max,self.ptdimAll,self.origsAll,self.csysesLocAll,self.DAll,self.SAll,self.RCAll,self.activeEdgeAll)
		print 'testlib.getUVW_at_points finshed'
		#return self.UVW
		iw=0
		vel=self.getFreeVelocity()
		for r in self.getRegionList():
			nPan=self.getRegion(r).getNumberOfPanels()
			ofsI=self.getPanelOffset(r)
			
			for i in range(0,nPan):	
				
				p=self.getRegion(r).getPanel(i)
				if p.isWake:
					if p.blWakeSolveVelocity:
						p.blUeUp=np.linalg.norm(self.UVW_Wake[iw]+self.getFreeVelocity())
						iw+=1
					else:
						p.blUeUp=np.linalg.norm(vel)
		#iw=self.nWPan

		for r in self.getRegionList():
			nPan=self.getRegion(r).getNumberOfPanels()
			ofsI=self.getPanelOffset(r)
			
			for i in range(0,nPan):	
				
				p=self.getRegion(r).getPanel(i)
				if p.isWake:
					if p.blWakeSolveVelocity:
						p.blUeDown=np.linalg.norm(self.UVW_Wake[iw]+self.getFreeVelocity())
						iw+=1
					else:
						p.blUeDown=np.linalg.norm(vel)
					
					
	def getVelocityAtPoints(self,pts):
		self.userPoi=pts.astype(np.float64)
		if self.userPoi.shape == (3,):
			self.userPoi=np.array([self.userPoi])
			npoi=1
		else:
			#self.userPoi.shape=(pts.shape[0]*3)
			npoi=pts.shape[0]

		self.userUVW=np.zeros((npoi,3),dtype=np.float64)
		npt=self.nPan
		ptdim_max=4

		self.allocPointsAll()
			
				
		#void getUVW_at_points(double*UVW, double*poi, int n, double*pts,int npt, int ptdim_max,int*ptdim, double*ORIGS, double*csysesLoc,double*D,double*S,double*RC){
		testlib.getUVW_at_points(self.userUVW,self.userPoi,npoi,self.ptsAll,self.nPan,ptdim_max,self.ptdimAll,self.origsAll,self.csysesLocAll,self.DAll,self.SAll,self.RCAll,self.activeEdgeAll)
		if pts.shape==(3,):
			return self.userUVW[0]+self.getFreeVelocity()
		else:
			return self.userUVW.reshape((pts.shape[0],3))+self.getFreeVelocity()


	def solveAij(self,normalsOrientation=1.0):
		self.buildGids()
		self.nWPan=self.nPan-self.nLPan
		print 'assembling Aij matrix'
		print strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
		if True:#self.nRuns < 1:
		#if self.nRuns < 1:
			self.B=np.zeros(self.nLPan,dtype=np.float64)
			self.BW=np.zeros(self.nLPan,dtype=np.float64)
			self.poi=np.zeros((self.nLPan,3),dtype=np.float64)

			self.pts=np.zeros((self.nLPan,3*4),dtype=np.float64)
			self.ptdim=np.ones(self.nLPan,dtype=np.int32)*4
			self.origs=np.zeros((self.nLPan,3),dtype=np.float64)
			self.csysesLoc=np.zeros((self.nLPan,9),dtype=np.float64)

			self.ptsW=np.zeros((self.nWPan,3*4),dtype=np.float64)
			self.ptdimW=np.ones(self.nWPan,dtype=np.int32)*4
			self.origsW=np.zeros((self.nWPan,3),dtype=np.float64)
			self.csysesLocW=np.zeros((self.nWPan,9),dtype=np.float64)
						
			
			
			self.D=np.ones(self.nLPan,dtype=np.float64)
			self.D0=np.zeros(self.nLPan,dtype=np.float64)
			self.DW=np.ones(self.nWPan,dtype=np.float64)
			self.SW=np.ones(self.nWPan,dtype=np.float64)
			
			self.S=np.ones(self.nLPan,dtype=np.float64)
			self.S0=np.zeros(self.nLPan,dtype=np.float64)
			self.A0=np.zeros(self.nLPan*self.nLPan ,dtype=np.float64)
			self.A=np.zeros((self.nLPan,self.nLPan) ,dtype=np.float64)
			self.gid=np.zeros(self.nPan,dtype=np.int32)


			self.firstPanel=np.ones(self.nWPan,dtype=np.int32)
			self.lastPanel=np.ones(self.nWPan,dtype=np.int32)

			self.isFirstPanel=np.zeros(self.nLPan,dtype=np.int32)
			self.isLastPanel=np.zeros(self.nLPan,dtype=np.int32)

			self.activeEdgeWake=np.zeros((self.nWPan,4),dtype=np.int32)


		else:
			self.B.fill(0.0)
			self.poi.fill(0.0)

			self.pts.fill(0.0)

			self.ptdim.fill(4)
			self.origs.fill(0.0)
			self.csysesLoc.fill(0.0)
			self.D.fill(1.0)
			self.DW.fill(1.0)
			self.D0.fill(0.0)
			self.S.fill(1.0)
			self.S0.fill(0.0)
			self.A.fill(0.0)
			self.A0.fill(0.0)
			self.gid.fill(0)


			self.ptsW.fill(0)
			self.ptdimW.fill(4)
			self.origsW.fill(0.0)
			self.csysesLocW.fill(0.0)


			self.firstPanel.fill(1)
			self.lastPanel.fill(1)
			self.isFirstPanel.fill(0)
			self.isLastPanel.fill(0)
			self.activeEdgeWake.fill(0)




		self.nRuns+=1
		n=self.nPan
		npt=self.nPan
		ptdim_max=4

		self.A0.shape = (self.nLPan*self.nLPan)
		print 'assembling RHS by testlib '+strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
		gidL=0
		for r in self.getRegionList():
			nWPan=self.getRegion(r).getNumberOfPanels()-self.getRegion(r).getNumberOfLiftPanels()
			nLPan=self.getRegion(r).getNumberOfLiftPanels()
			ofsI=self.getPanelOffset(r)
			ofsIL=self.getLiftPanelOffset(r)
			
			for i in range(0,nLPan):
				p=self.getRegion(r).getLiftPanel(i)
				self.poi[i+ofsIL]=np.float64(p.controlPoint)
				# OFF BODY KUTTA CONDITION
				#if False:#self.blIter>0:
					#if p.isFirst or p.isLast:
						#self.poi[i+ofsIL]=np.float64(p.getCpoint()- p.getNV()*p.displThick)

				if p.isFirst:
					self.isFirstPanel[i+ofsIL]=1
				if p.isLast:
					self.isLastPanel[i+ofsIL]=1
					
			for i in range(0,nLPan):
				p=self.getRegion(r).getLiftPanel(i)
				p.gidL=gidL
				gidL+=1
				self.pts[i+ofsIL]=np.float64( p.pts.flatten())
				self.origs[i+ofsIL]=np.float64(p.controlPoint)
				self.csysesLoc[i+ofsIL]=np.float64(p.csys.flatten())

			for i in range(0,nLPan):	
				p=self.getRegion(r).getLiftPanel(i)
				if p.isWake:
					self.S[i+ofsIL]=np.float64(0.0-p.getBlowingVelocity())
				else:
					rvector=p.getCpoint()-self.COR
					vrot=np.cross(self.omega,rvector)
					self.S[i+ofsIL]=np.float64(np.dot((self.getFreeVelocity()+vrot+p.vFarWake+p.vRel) ,p.normalVect)-p.getBlowingVelocity())

		iw=0
		for r in self.getRegionList():
			nPan=self.getRegion(r).getNumberOfPanels()
			ofsI=self.getPanelOffset(r)
			for i in range(0,nPan):	
				p=self.getRegion(r).getPanel(i)
				if p.isWake:
					#if p.G==0.0:
						#p.activeEdges*=0
					# trailing edge potential jump, Youngren 1983 - QUADPAN
					pFirst=self.getPanelByGid(p.firstPanel)
					pLast=self.getPanelByGid(p.lastPanel)
					rte=pFirst.getCpoint()-pLast.getCpoint()
					if self.blIter>0:
						rte=(pFirst.getBlCpoint(1.)- pLast.getBlCpoint(1.))
						rte=(pFirst.getBlCpoint(0.)- pLast.getBlCpoint(0.))

					bv=0.5*(pFirst.dragVect + pLast.dragVect)
					bv/=np.linalg.norm(bv)

					xn=0.5*(pFirst.csys[1] + pLast.csys[1])
					xn/=np.linalg.norm(xn)
					
					nb=np.cross(bv,xn)
						
					self.DW[iw]=np.dot(self.getFreeVelocity(),nb)*np.dot(rte,nb)
					#self.SW[iw]=-p.w
					
					SFirst=np.float64(np.dot((self.getFreeVelocity()+vrot) ,pFirst.normalVect)-pFirst.getBlowingVelocity())
					SLast=np.float64(np.dot((self.getFreeVelocity()+vrot) ,pLast.normalVect)-pLast.getBlowingVelocity())
					
					p.S=-p.getBlowingVelocity()
					self.SW[iw]=p.getBlowingVelocity()
					
					self.ptsW[iw]=np.float64( p.pts.flatten())
					self.origsW[iw]=np.float64(p.controlPoint)
					self.csysesLocW[iw]=np.float64(p.csys.flatten())
					pFirst=self.getPanelByGid(p.firstPanel)
					pLast=self.getPanelByGid(p.lastPanel)
					self.firstPanel[iw]=pFirst.gidL
					self.lastPanel[iw]=pLast.gidL
					self.activeEdgeWake[iw]=p.activeEdges
					iw+=1
			
			

		affinity.set_process_affinity_mask(0,2**multiprocessing.cpu_count()-1)
		testlib.getPHI_at_points(self.B,self.poi,self.nLPan,self.pts,self.nLPan,ptdim_max,self.ptdim,self.origs,self.csysesLoc,self.D0,self.S)

		testlib.getPHI_at_points(self.BW,self.poi,self.nLPan,self.ptsW,self.nWPan,ptdim_max,self.ptdim,self.origsW,self.csysesLocW,self.DW,self.SW)
		

		print 'assembling AIJ by testlib '+strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
		affinity.set_process_affinity_mask(0,2**multiprocessing.cpu_count()-1)
		testlib.getAij_matrix1(self.A0,self.poi,self.nLPan,self.pts,self.origs,self.csysesLoc,self.nWPan,self.ptsW,self.origsW,self.csysesLocW,self.firstPanel,self.lastPanel,self.isFirstPanel,self.isLastPanel,self.activeEdgeWake)
		

		self.A0.shape = (self.nLPan,self.nLPan)

		##self.B*=-1.0#
		self.B=-self.B
		self.B-=self.BW
		
		
		print 'Aij matrix assembled'
		print 'solving Aij'
		affinity.set_process_affinity_mask(0,2**multiprocessing.cpu_count()-1)
		self.G=lgmres(self.A0,self.B,tol=1e-9)[0]
		print 'postprocessing results'
		print strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())

		for r in self.getRegionList():
			nLPan=self.getRegion(r).getNumberOfLiftPanels()
			ofs=self.getLiftPanelOffset(r)
			region=self.getRegion(r)
			for i in range(0,nLPan):
				region.getLiftPanel(i).setG(self.G[ofs+i])
				region.getLiftPanel(i).setS(self.S[ofs+i])


		for i in self.gids:
			p=self.getPanelByGid(i)
			if p.isFirst:
				g=p.G- self.getPanelByGid(p.lastPanel).G
				for j in p.wakes:
					self.getPanelByGid(j).G=np.copy(g) ## !!! Add (-) just to adjust velocity survey, check trailing edge singularity
			#if np.all(p.activeEdges==0):
					#p.G=0.
					
					
		for r in self.getRegionList():
			nPan=self.getRegion(r).getNumberOfPanels()
			region=self.getRegion(r)
			for i in range(0,nPan):	
				region.getPanel(i).setForce()
				
		dx=['B','BW','poi','pts','ptdim','origs','csysesLoc','ptsW','ptdimW','origsW','csysesLocW','D','D0','DW','SW','S','S0','A0','A','gid','firstPanel','lastPanel','isFirstPanel','isLastPanel','activeEdgeWake']
		for i in dx:
			delattr(self,i)
			setattr(self,i,None)
		ctypes._reset_cache()
		print 'results computed: '+strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())

	def relaxWake(self,w,restricted,KW=0.05):
		pass
		#ptsW=np.zeros((self.nWPan*4,3))
		#displW=np.zeros((self.nWPan*4,3))
		#iw=0
		#ir=0
		#for i in range(0,self.nPan):
			#p=self.getPanelByGid(i)
			#if p.isWake:
				#for j in range(0,4):
					#ptsW[iw*4+j]=p.pts[j]
					#displW[iw*4+j]=self.getVelocityAtPoint(ptsW[iw*4+j])-np.array([self.freeVel[0],0,0])
				#iw+=1
			#if p.isFirst:
				#for j in range(0,4):
					#restricted[ir*4+j]=p.pts[j]
				#ir+=1

		#for iw in range(0,ptsW.shape[0]):
			#for ir in range(0,restricted.shape[0]):
				#if np.linalg.norm(restricted[ir]-ptsW[iw]) <0.0001:
					#displW[iw]*=0.

					
		#iw=0
		#ir=0
		#for i in range(0,self.nPan):
			#p=self.getPanelByGid(i)
			#if p.isWake:
				#for j in range(0,4):
					#p.pts[j]+=displW[iw*4+j]*KW
				#iw+=1
	
	def setWakePanelsGByBoundingBoxes(self,bbs,G=0.0):
		for bb in bbs:
			gidW=self.getPanelGidsByBoundingBox(bb[0],bb[1])
			for j in gidW:
				p=self.getPanelByGid(j)
				if p.isWake:
					p.activeEdges*=0
					#p.G=G

class VlmProblem(GlobalEnvironment):
	
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


	def getSolution0(self):
		dom1=self.dom1
		dom1.setPanelsFreeStream(self.getFreeVelocity())
		dom1.reinitPanelsGeometry()
		#dom1.setPanelsG(1.0)
		dom1.solveAij0()

	def getSolutionAeroelastic(self,stopF,stopDim=None,fullRun=False,err=1e-3,maxIt=10,reset=True,wake=None,dt=0.0):
		ii=np.linspace(0,maxIt-1,maxIt)
		fi=np.zeros(maxIt)
		if self.fem.lsCount==0:
			self.getSolution()
		for i in range(0,maxIt):
			print('solving aeroelastic iteration '+str(i))

			if stopDim is None:
				stop=stopF()
			else:
				stop=stopF()[stopDim]
			stopOld=copy.copy(stop)
			
			if dt!=0.0:
				self.fem.computeKinematics(dt)
				self.fem.setDisplacementHistory()
			try:
				self.elastic.mapInertias(self.fem,self.inertia)
			except:
				print('Skipping inertia forces')
				pass
			force0=copy.copy(self.elastic.force2Fem())
			if reset and i ==0:
				print('restarting aeroelastic')
				self.fem.force=force0
				self.fem.processInputs()
				self.fem.solveDisplacementsNonlinear(nls=3,tol=1e-3,maxIt=20)
			else:
				self.fem.solveLoadStep(self.fem.forceDict2Vector(force0),tol=1e-3)
			
			self.elastic.deform(self.fem.DISPLACEMENT_LS.reshape(self.fem.NODE.shape))

			#self.elastic.buildCsys()
			self.getSolution()
			if wake is not None:
				wake.relax(self.dom1,0.5)
			if stopDim is None:
				stop=stopF()
			else:
				stop=stopF()[stopDim]
			fi[i]=stop
			print 'stop: '+str(stop)+' stopOld: '+str(stopOld)
			if np.abs(stop - stopOld) < err and not fullRun:
				break
			#self.elastic.buildCsys()
		return ii,fi
#			print('Force: '+str(wing.getForceAtRefPT()))
#	print fem.F.toarray()[6:42]
#			f0=wing.getForceAtRefPT()
#			f=f0#/0.5/rho/sp**2.0/area
#			cm=f[5]
#			cl=f[1]
#			cd=f[0]
#			cl=-f[0]*np.sin(alpha)+f[1]*np.cos(alpha)
#			cd= f[0]*np.cos(alpha)+f[1]*np.sin(alpha)
#			cm=f[5]
#			fi[i]=fem.D[0][1]

	def getSolutionWithFarWake(self,farWake=None,useLastStep=False,vortexRC=0.1,relaxingWake=None,kw=0.1,nWakeRelax=1,wakeFreeDOF=np.array([1,1,1])):
		if farWake is not None:
			wakeFar=GlobalDomain()
			wakeFar.addRegion(farWake)
			wakeFar.vortexRC=vortexRC
			for i in wakeFar.gids:
				p=wakeFar.getPanelByGid(i)
				p.activeEdges.fill(1)
		if not(self.dom1.nRuns>0  and useLastStep):		
			self.getSolution()
			

		if farWake is not None:
			poiHT={}
			j=0
			for i in self.dom1.gids:
				p=self.dom1.getPanelByGid(i)
				if not p.isWake:
					poiHT[j]=i
					j+=1
			n=j
			poi=np.zeros((n,3))
			for i in range(0,n):
				poi[i]=self.dom1.getPanelByGid(poiHT[i]).getCpoint()

			vff=wakeFar.getVelocityAtPoints(poi)-self.getFreeVelocity()
			for i in range(0,n):
				p=self.dom1.getPanelByGid(poiHT[i])
				p.vFarWake=vff[i]
			
			
			for i in wakeFar.gids:
				p=self.dom1.getPanelByGid(i)
				p.activeEdges.fill(0)

			self.getSolution()

			if relaxingWake is not None:
				#for k in relaxingWake:
					#relaxingWake[k].setHorseshoeWake(self.dom1)
				for j in range(0,nWakeRelax):
					for k in relaxingWake:
						relaxingWake[k].relax(self.dom1, kw=kw, findDiplacement=True, deform=False)
						
					for k in relaxingWake:
						relaxingWake[k].relax(self.dom1, kw=kw, findDiplacement=False, deform=True,freeDOF=wakeFreeDOF)

			for i in wakeFar.gids:
				p=self.dom1.getPanelByGid(i)
				p.activeEdges.fill(1)
			
			
	def resetBLBlowing(self):
		for r in self.dom1.getRegionList():
			for i in self.dom1.regions[r].getPanelList():
				#self.dom1.regions[r].getPanel(i).blIter=0
				self.dom1.regions[r].blIter=0
				self.dom1.blIter=0
				self.blIter=0
				GlobalEnvironment.blIter=0
				self.dom1.regions[r].getPanel(i).Ue=0.0
				self.dom1.regions[r].getPanel(i).w=0.0
				self.dom1.regions[r].getPanel(i).vnpbl=0.0
				self.dom1.regions[r].getPanel(i).dlspbl=0.0
				self.dom1.regions[r].getPanel(i).ue=0.0
				self.dom1.regions[r].getPanel(i).uebl=0.0
				self.dom1.regions[r].getPanel(i).uedispbl=0.0


				self.dom1.regions[r].getPanel(i).blDtUp=0.
				self.dom1.regions[r].getPanel(i).blDtDown=0.
				self.dom1.regions[r].getPanel(i).blUeUp=0.
				self.dom1.regions[r].getPanel(i).blUeDown=0.
				self.dom1.regions[r].getPanel(i).blVnUp=0.0
				self.dom1.regions[r].getPanel(i).blVnDown=0.
	
	def getSolutionWithBL(self,pgs,maxIt=30,maxSweeps=10,stopFunc=None,sfArgs={},tol=2.):
		fcs=np.zeros(len(pgs))
		fcsNew=np.zeros(len(pgs))

		#for r in self.dom1.getRegionList():
			#for k in self.dom1.regions[r].getPanelList():
				#self.dom1.regions[r].getPanel(k).w=0.0
				#self.dom1.regions[r].getPanel(k).vnpbl=0.0
		self.getSolution()
		
		for i in range(0,maxIt):
				
			j=0
			self.dom1.setVelocityAtWakeBlPoints()
			#for r in self.dom1.getRegionList():
				#for k in self.dom1.regions[r].getPanelList():
					#self.dom1.regions[r].getPanel(k).w=0.0
					#self.dom1.regions[r].getPanel(k).vnpbl=0.0

			for pg in pgs:
				for nsw in range(0,maxSweeps):
					pg.examineBLAnalitic(self.dom1)
					pg.blIter+=1
			for r in self.dom1.getRegionList():
				for k in self.dom1.regions[r].getPanelList():
					pan=self.dom1.regions[r].getPanel(k)
					if pan.isWake:
						pan.w=(-pan.blVnUp + pan.blVnDown)*-1.
			self.getSolution()
			fcsNew[j]=stopFunc(sfArgs)
			print fcsNew
			if np.all( np.abs(fcsNew-fcs) < tol) and i >2:
				break
			fcs=np.copy(fcsNew)
			self.dom1.blIter+=1
			GlobalEnvironment.blIter+=1

			#for pg in pgs:
				#pg.setBlowing()
			
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
				if not panelI.isWake or plotWake:
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
		try:
			for i in self.elastic.interactors:
				I=self.elastic.interactors[i]
				self.view.addCSYS(I.orig,I.csysLoc,scaling=0.2)
		except:
			pass
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
	
	
		