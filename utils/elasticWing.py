import numpy as np
import copy
from utils.mesh import *
from panelMethod.panelMethod1 import *
from panelMethod.env import GlobalEnvironment
from utils.b3Vect import B3vect 
from utils.bCad import BCad

from fem.fem4 import FEM1
from fem.beam3d_nl import Beam3D
from utils.interactor import PanelBeamInteractor
import utils.loftedSurface as ls
import utils.csys as cs

from utils.userFunctions import *

class ApparrentMass(object):
	def __init__(self):
		self.acc=np.zeros(3)
		self.p=np.zeros(3)
		self.m=0.0
	
	def getForceAtPT(self,pt):
		acc=self.acc
		f=np.zeros(6)
		f[0:3]=-acc*self.m
		f[3:6]=np.cross(self.p-pt, -acc*self.m)
		return f
		
class ElasticWing1(object):
	def __init__(self):
		self.refPT=np.zeros(3)
		self.cuts={}
		self.interactors={}
		self.nSubSec=0
		self.pts={}
		pass

	def setRefPT(self,PT):
		self.refPT=PT

	def getSubSectionSparNodes(self,section,spar,j1=0,j2=2,oi=0,i1=20,i2=40):
		nsec=section.nSpanwise
		npan=section.nPanChordwise
		pts1=findSubsectionSparPoints(section,spar,j1,oi,i1,i2)
		pts2=findSubsectionSparPoints(section,spar,j2,oi,i1,i2)
		return np.vstack((pts2,pts1[-1,:]))


	def makeInteractors1(self,section,nodes,ofsI=0,reverse=False,integration='12',swapWake=False,oVector=np.array([0,0,1.])):
		nsec=section.nSpanwise
		npan=section.nPanChordwise
		self.nSubSec+=1
		pgs=extractSpanwisePanelGroupsFromSection(section)
		kk=pgs.keys()
		if reverse:
			kk=kk[::-1]
		for i in range(0,len(kk)):
			pg=pgs[kk[i]]
			j=len(self.interactors)
			self.interactors[j]=PanelBeamInteractor()
			self.interactors[j].integration=integration
			self.interactors[j].section[0]=ls.LoftSection()
			self.interactors[j].section[1]=ls.LoftSection()
			self.interactors[j].section[0].orig=nodes[i].copy()
			self.interactors[j].section[1].orig=nodes[i+1].copy()
			self.interactors[j].oVector=oVector

			if integration=='21':
				self.interactors[j].node[0]=nodes[i]
				self.interactors[j].node[1]=nodes[i+1]
				self.interactors[j].ids[0]=i+ofsI
				self.interactors[j].ids[1]=i+1+ofsI
			if integration=='12':
				self.interactors[j].node[0]=nodes[i]
				self.interactors[j].node[1]=nodes[i+1]
				self.interactors[j].ids[0]=i+ofsI
				self.interactors[j].ids[1]=i+1+ofsI
			
			for s in pg.getPanelList():
				if integration=='12':
					c1=self.interactors[j].section[1]
					c2=self.interactors[j].section[0]
				if integration=='21':
					c1=self.interactors[j].section[0]
					c2=self.interactors[j].section[1]
					
				pan=pg.getPanel(s)
				#c1.addPoint(pan.pts[0])
				#c1.addPoint(pan.pts[3])
				#c2.addPoint(pan.pts[1])
				#c2.addPoint(pan.pts[2])
				
				if not pan.isWake:
				#c1.addSurfacePoint(pan.pts[0])
				#c1.addSurfacePoint(pan.pts[3])
				#c2.addSurfacePoint(pan.pts[1])
				#c2.addSurfacePoint(pan.pts[2])

					c1.addPoint(pan.pts[0])
					c1.addPoint(pan.pts[3])
					c2.addPoint(pan.pts[1])
					c2.addPoint(pan.pts[2])
						
				if pan.isWake and swapWake:
					c1.addPoint(pan.pts[1])
					c1.addPoint(pan.pts[2])
					c2.addPoint(pan.pts[0])
					c2.addPoint(pan.pts[3])
				
				if pan.isWake and not swapWake:
					c1.addPoint(pan.pts[0])
					c1.addPoint(pan.pts[3])
					c2.addPoint(pan.pts[1])
					c2.addPoint(pan.pts[2])
			
			self.interactors[j].setPanelGroup(pg)

			if integration == '12':
				self.interactors[j].pg.refPT=nodes[i+1].copy()
				#self.interactors[j].orig=nodes[i+1].copy()
				self.interactors[j].csysPT1=nodes[i+1]
			if integration == '21':
				self.interactors[j].pg.refPT=nodes[i].copy()
				#self.interactors[j].orig=nodes[i].copy()
				self.interactors[j].csysPT1=nodes[i]
			self.interactors[j].orig= 0.5*(self.interactors[j].node[0]+self.interactors[j].node[1])


	def addSubSection(self,section,nodes,spar,ofsCC=0,integration='12',j1=0,j2=2,oi=0,i1=20,i2=40,oVector=np.array([0,0,1.])):
		
		nsec=section.nSpanwise
		npan=section.nPanChordwise
		pts1=findSubsectionSparPoints(section,spar,j1,oi,i1,i2)
		pts2=findSubsectionSparPoints(section,spar,j2,oi,i1,i2)
		self.pts[self.nSubSec]={}
		self.pts[self.nSubSec][1]=pts1
		self.pts[self.nSubSec][2]=pts2
		self.nSubSec+=1
		ofsC=len(self.cuts)+ofsCC
		ofsI=len(self.interactors)
		pgs=extractSpanwisePanelGroupsFromSection(section)
		for i in pgs:
			self.interactors[i+ofsI]=PanelBeamInteractor()
			self.interactors[i+ofsI].integration=integration
			self.interactors[i+ofsI].section[0]=ls.LoftSection()
			self.interactors[i+ofsI].section[1]=ls.LoftSection()
			if integration == '12':
				self.interactors[i+ofsI].orig=pts1[i]
				self.interactors[i+ofsI].section[0].orig=pts1[i]
				self.interactors[i+ofsI].section[1].orig=pts2[i]
			if integration == '21':
				self.interactors[i+ofsI].orig=pts2[i]
				self.interactors[i+ofsI].section[0].orig=pts1[i]
				self.interactors[i+ofsI].section[1].orig=pts2[i]
			
			for s in pgs[i].getPanelList():

				c1=self.interactors[i+ofsI].section[0]
				c2=self.interactors[i+ofsI].section[1]
				pan=pgs[i].getPanel(s)
				c2.addPoint(pan.pts[0])
				c2.addPoint(pan.pts[3])
				c1.addPoint(pan.pts[1])
				c1.addPoint(pan.pts[2])
				
				if not pan.isWake:
					c2.addSurfacePoint(pan.P1)
					c1.addSurfacePoint(pan.P2)
				if pan.isLast:
					c2.addSurfacePoint(pan.P4)
					c1.addSurfacePoint(pan.P3)
			

			if not pan.isWake:
				c2.addSurfacePoint(pan.P2)
				c1.addSurfacePoint(pan.P3)
			
			self.interactors[i+ofsI].setPanelGroup(pgs[i])


			#if integration == '12':
				#self.interactors[i+ofsI].integration = '12'
				#self.interactors[i+ofsI].orig=self.interactors[i+ofsI].section[1].orig
			#if integration == '21':
				#self.interactors[i+ofsI].orig=self.interactors[i+ofsI].section[0].orig
				#self.interactors[i+ofsI].integration = '21'
			self.interactors[i+ofsI].setOrientationVector(oVector)
		#return pgs,cuts

	def connectNodes(self,nodes,start,stop,ofsN=0): 
		interactors=self.interactors
		for i in range(start,stop):
						

			nodes[i+ofsN]=interactors[i].section[0].orig.copy()
			nodes[i+ofsN+1]=interactors[i].section[1].orig.copy()

			#interactors[i].section[0].orig=nodes[i+ofsN]
			#interactors[i].section[1].orig=	nodes[i+ofsN+1]
			interactors[i].node[0]=nodes[i+ofsN]
			interactors[i].node[1]=nodes[i+ofsN+1]
			interactors[i].ids[0]=i+ofsN
			interactors[i].ids[1]=i+1+ofsN
			##if i >0:
				##interactors[i-1].section[1].csys=interactors[i].section[0].csys
				##interactors[i-1].section[1].orig=nodes[i+ofsN]
				##interactors[i-1].node[1]=nodes[i+ofsN]
				##interactors[i-1].ids[1]=i+ofsN

			if interactors[i].integration == '12':
				interactors[i].orig=interactors[i].node[1].copy()
				interactors[i].csysPT1=interactors[i].node[0]
			if interactors[i].integration == '21':
				interactors[i].orig=interactors[i].node[0]
				interactors[i].csysPT1=interactors[i].node[1].copy()

			if interactors[i].integration == '12':
				interactors[i].pg.refPT=nodes[i+ofsN+1].copy()
			if interactors[i].integration == '21':
				interactors[i].pg.refPT=nodes[i+ofsN].copy()

			
			self.initialNodes=copy.copy(nodes)
			self.nodes=nodes
			
			
		
	def setBySections(self,section,nodes,rootCH025,ofsY,symCoor=2):
		ns=len(section)
		cuts=self.cuts
		interactors=self.interactors

		for i in section:
			interactors[i]=PanelBeamInteractor()
			cuts[i*2]=ls.LoftSection()
			cuts[i*2+1]=ls.LoftSection()
			for s in section[i].getPanelList():
				pan=section[i].getPanel(s)
				cuts[i*2].orig=np.array([rootCH025,0,section[i].panels[s].pts[0][2]])
				cuts[i*2].addPoint(pan.P1)
				

				cuts[i*2+1].orig=np.array([rootCH025,0,section[i].panels[s].pts[1][2]])
				cuts[i*2+1].addPoint(pan.P2)
				if not pan.isWake:
					cuts[i*2].addSurfacePoint(pan.P1)
					cuts[i*2+1].addSurfacePoint(pan.P2)
			cuts[i*2].addPoint(pan.P4)
			cuts[i*2+1].addPoint(pan.P3)
			if not pan.isWake:
				cuts[i*2].addSurfacePoint(pan.P4)
				cuts[i*2+1].addSurfacePoint(pan.P3)

			interactors[i].setPanelGroup(section[i])
			interactors[i].setSection(0,cuts[i*2]);
			interactors[i].setSection(1,cuts[i*2+1])


		for i in interactors:
			if interactors[i].section[1].orig[symCoor] > 0.0:
				interactors[i].integration = '21'


			if interactors[i].integration == '12':
				interactors[i].orig=interactors[i].section[1].orig
			if interactors[i].integration == '21':
				interactors[i].orig=interactors[i].section[0].orig
			
			interactors[i].setOrientationVector(ofsY)

			nodes[i]=interactors[i].section[0].orig
			nodes[i+1]=interactors[i].section[1].orig

			interactors[i].section[0].orig=nodes[i]
			interactors[i].section[1].orig=	nodes[i+1]
			interactors[i].node[0]=nodes[i]
			interactors[i].node[1]=nodes[i+1]
			interactors[i].ids[0]=i
			interactors[i].ids[1]=i+1
			if i >0:
				interactors[i-1].section[1].csys=interactors[i].section[0].csys
				interactors[i-1].section[1].orig=nodes[i]
				interactors[i-1].node[1]=nodes[i]
				interactors[i-1].ids[1]=i

			if interactors[i].integration == '12':
				interactors[i].orig=interactors[i].node[1]
				interactors[i].csysPT1=interactors[i].node[0]
			if interactors[i].integration == '21':
				interactors[i].orig=interactors[i].node[0]
				interactors[i].csysPT1=interactors[i].node[1]

			if interactors[i].integration == '12':
				interactors[i].pg.refPT=nodes[i+1]
			if interactors[i].integration == '21':
				interactors[i].pg.refPT=nodes[i]

			if interactors[i].section[1].orig[symCoor] <= 0.0:
				for j in range(0,i):
					interactors[i].addPreviousInteractor(interactors[j])

			if interactors[i].section[1].orig[symCoor] > 0.0:
				for j in range(ns-1,i,-1):
					interactors[i].addPreviousInteractor(interactors[j])
			self.initialNodes=copy.copy(nodes)
			self.nodes=nodes


		
	def setBySectionsWithSpine(self,section,nodes,spine,ofsY,symCoor=2,symAxisFlip=False):
		ns=len(section)
		cuts=self.cuts
		interactors=self.interactors
		spl=copy.copy(spine)
		for i in section:
			interactors[i]=PanelBeamInteractor()
			cuts[i*2]=ls.LoftSection()
			cuts[i*2+1]=ls.LoftSection()

			l=section[i].panels[0].pts[0][symCoor]
			xi=np.interp(l,spl[:,symCoor],spl[:,0])
			yi=np.interp(l,spl[:,symCoor],spl[:,1])
			zi=np.interp(l,spl[:,symCoor],spl[:,2])
			cuts[i*2].orig=np.array([xi,yi,zi])
			cuts[i*2].orig[symCoor]=l


			l=section[i].panels[0].pts[1][symCoor]
			xi=np.interp(l,spl[:,symCoor],spl[:,0])
			yi=np.interp(l,spl[:,symCoor],spl[:,1])
			zi=np.interp(l,spl[:,symCoor],spl[:,2])
			cuts[i*2+1].orig=np.array([xi,yi,zi])
			cuts[i*2+1].orig[symCoor]=l


			for s in section[i].getPanelList():
				pan=section[i].getPanel(s)
				
				cuts[i*2].addPoint(pan.P1)
				

				cuts[i*2+1].addPoint(pan.P2)
				if not pan.isWake:
					cuts[i*2].addSurfacePoint(pan.P1)
					cuts[i*2+1].addSurfacePoint(pan.P2)
			cuts[i*2].addPoint(pan.P4)
			cuts[i*2+1].addPoint(pan.P3)
			if not pan.isWake:
				cuts[i*2].addSurfacePoint(pan.P4)
				cuts[i*2+1].addSurfacePoint(pan.P3)

			interactors[i].setPanelGroup(section[i])
			interactors[i].setSection(0,cuts[i*2]);
			interactors[i].setSection(1,cuts[i*2+1])


		for i in interactors:
			if interactors[i].section[1].orig[symCoor] > 0.0 and (not symAxisFlip):
				interactors[i].integration = '21'


			if interactors[i].integration == '12':
				interactors[i].orig=interactors[i].section[1].orig
			if interactors[i].integration == '21':
				interactors[i].orig=interactors[i].section[0].orig
			
			interactors[i].setOrientationVector(ofsY)

			nodes[i]=interactors[i].section[0].orig
			nodes[i+1]=interactors[i].section[1].orig

			interactors[i].section[0].orig=nodes[i]
			interactors[i].section[1].orig=	nodes[i+1]
			interactors[i].node[0]=nodes[i]
			interactors[i].node[1]=nodes[i+1]
			interactors[i].ids[0]=i
			interactors[i].ids[1]=i+1
			if i >0:
				interactors[i-1].section[1].csys=interactors[i].section[0].csys
				interactors[i-1].section[1].orig=nodes[i]
				interactors[i-1].node[1]=nodes[i]
				interactors[i-1].ids[1]=i

			if interactors[i].integration == '12':
				interactors[i].orig=interactors[i].node[1]
				interactors[i].csysPT1=interactors[i].node[0]
			if interactors[i].integration == '21':
				interactors[i].orig=interactors[i].node[0]
				interactors[i].csysPT1=interactors[i].node[1]

			if interactors[i].integration == '12':
				interactors[i].pg.refPT=nodes[i+1]
			if interactors[i].integration == '21':
				interactors[i].pg.refPT=nodes[i]

			if not symAxisFlip:
				if interactors[i].section[1].orig[symCoor] <= 0.0:
					for j in range(0,i):
						interactors[i].addPreviousInteractor(interactors[j])
	
				if interactors[i].section[1].orig[symCoor] > 0.0:
					for j in range(ns-1,i,-1):
						interactors[i].addPreviousInteractor(interactors[j])
	
			if symAxisFlip:
				if interactors[i].section[1].orig[symCoor] >= 0.0:
					for j in range(0,i):
						interactors[i].addPreviousInteractor(interactors[j])

#				if interactors[i].section[1].orig[symCoor] <= 0.0:
#					for j in range(ns-1,i,-1):
#						interactors[i].addPreviousInteractor(interactors[j])


			self.initialNodes=copy.copy(nodes)
			self.nodes=nodes


	def setRefDistances(self,PT):
		interactors=self.interactors
		for i in interactors:
			interactors[i].refDist=np.linalg.norm( interactors[i].orig - PT )

	def force2Fem(self):
		interactors=self.interactors
		ff0={}
		fn=[]
		for i in interactors:
			ff0[i]=interactors[i].forces2fem()
			fn.append(ff0[i]['node'])
		f={}
		addedN=[]
		for i in ff0:
			if not fn[i] in addedN:
				f[i]=ff0[i]
			else:
				for j in f:
					if f[j]['node'] == fn[i]:
						f[j]['force']+=ff0[i]['force']
			addedN.append(fn[i])
		return f


	def deform(self,D,moveOPoints=True):
		interactors=self.interactors
		for i in interactors:
			interactors[i].fem2deform(D)
		#self.buildCsys()

	def _rotateOPoint(self,nodes,D):
		interactors=self.interactors
		for i in interactors:
			interactors[i]._rotateOPoint(nodes,D)


	def _rotateOVector(self,nodes,D):
		interactors=self.interactors
		for i in interactors:
			interactors[i]._rotateOVector(nodes,D)
	

	def buildCsys(self):
		interactors=self.interactors
		for i in interactors:
			interactors[i].buildCsys()


	def getLocalCumulativeForces(self,j):
		interactors=self.interactors
		x=np.zeros(len(interactors)); y=np.zeros(len(interactors));
		for i in interactors:
			I=interactors[i]; 
			x[i]=I.refDist;
			if I.integration == '12':
				x[i]=-I.refDist; 
			f=cs.transformVector(I.getCumulativeForceAtRefPT()[0:3], I.csysLoc, I.csys)
			m=cs.transformVector(I.getCumulativeForceAtRefPT()[3:6], I.csysLoc, I.csys)
			F=np.hstack((f,m))
			y[i]= F[j];
		return x,y

	def getCumulativeForces(self,j):
		interactors=self.interactors
		x=np.zeros(len(interactors)); y=np.zeros(len(interactors));
		for i in interactors:
			I=interactors[i]; 
			x[i]=I.refDist;
			if I.integration == '12':
				x[i]=-I.refDist; 
			f=I.getCumulativeForceAtRefPT()[0:3]
			m=I.getCumulativeForceAtRefPT()[3:6]
			F=np.hstack((f,m))
			y[i]= F[j];
		return x,y

	def getLocalForces(self,j):
		interactors=self.interactors
		x=np.zeros(len(interactors)); y=np.zeros(len(interactors));
		for i in interactors:
			I=interactors[i]; 
			x[i]=I.refDist;
			if I.integration == '12':
				x[i]=-I.refDist;
			f=cs.transformVector(I.getForceAtRefPT()[0:3], I.csysLoc, I.csys)
			m=cs.transformVector(I.getForceAtRefPT()[3:6], I.csysLoc, I.csys)
			F=np.hstack((f,m))
			y[i]= F[j];
		return x,y




	def resetGeometry(self):
		for i in range(0,self.nodes.shape[0]):
			for j in range(0,self.nodes.shape[1]):
				self.nodes[i][j]=self.initialNodes[i][j]
		for i in self.interactors:
			self.interactors[i].resetGeometry()
			self.interactors[i].buildCsys()
			for j in self.interactors[i].pg.getPanelList():
				self.interactors[i].pg.getPanel(j).initGeometry()
		
	def exportGeometry(self,folder,limit=None):
		
		g={}
		j=0
		for i in self.interactors:
			if limit == None:
				g[j]= self.interactors[i].section[0].pts2arr()
			else:
				g[j]= self.interactors[i].section[0].pts2arr()[0:limit,:]
			j+=1
		if limit == None:
			g[j]=self.interactors[i].section[1].pts2arr()
		else:
			g[j]= self.interactors[i].section[1].pts2arr()[0:limit,:]

		for i in g:          
			np.savetxt(folder+'/section'+str(i)+'.dat',g[i])

	def mapInertias(self,fem,im):
		acc_NODE=fem.ACCELERATION.reshape(fem.NODE.shape)
		for i in self.interactors:
			ia=self.interactors[i]
			if ia.integration == '12':
				c=ia.ids[1]
			if ia.integration=='21':
				c=ia.ids[0]
			accN=acc_NODE[c][0:3]
			
			for m in ia.mass:
				m.acc=im.getAccelerationAtPT(m.p)+accN
	
	def makeInertiaPanels(self,rho,thick):
		for i in self.interactors:
			ia=self.interactors[i]
			for j in ia.pg.panels:
				p=ia.pg.panels[j]
				if not p.isWake:
					mm=ApparrentMass()
					mm.p=np.copy(p.getCpoint())
					mm.m=p.getArea()*rho*thick
					ia.mass.append(copy.copy(mm))
					
	def makeInertiaSpar(self,rho):
		for i in self.interactors:
			ia=self.interactors[i]
			mm=ApparrentMass()
			mm.p=ia.getCentroid()
			mm.m=ia.getNodeDistace()*rho
			ia.mass.append(copy.copy(mm))
			
	def getMass(self):
		mass=0.
		for i in self.interactors:
			ia=self.interactors[i]
			for m in ia.mass:
				mass+=m.m
		return mass
	
	def getCG(self):
		mass=0.
		moment=np.zeros(3)
		for i in self.interactors:
			ia=self.interactors[i]
			for m in ia.mass:
				mass+=m.m
				moment+=m.p*m.m
		return moment/mass
	