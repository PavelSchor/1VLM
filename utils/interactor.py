import numpy as np
import copy
from utils.bCad import BCad
import utils.csys as cs
bc=BCad()


class PanelBeamInteractor(object):
	def __init__(self):
		self.orientation='byVector'
		self.csys=np.array([ [1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
		self.csysLoc=np.array([ [1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
#		self.orig=np.zeros(3)
#		self.oPoint=np.zeros(3)
		self.previousPanelGroups=[]
		self.integration='12'
		self.ids={}
		self.codes={}
		self.section={}
		self.force={}
		self.deformation={}
		self.node={}
		self.oVector=np.array([0.,1.,0.])
		self.oVectorInitial=np.array([0.,1.,0.])
		self.mass=[]
		pass

	def getCentroid(self):
		return (self.getNode(0)+self.getNode(1))*0.5

	def getNodeDistace(self):
		return np.linalg.norm(self.getNode(0)-self.getNode(1))
	
	def setOrientationVector(self,v):
		for i in range(0,3):
			self.oVector[i]=v[i]
			self.oVectorInitial[i]=v[i]

	def addPreviousInteractor(self,pg):
		self.previousPanelGroups.append(pg)

	def getForceAtPT(self,PT):
		return self.pg.getForceAtPT(PT)

	def getForceAtRefPT(self):
		return self.pg.getForceAtPT(self.pg.refPT)

	def getCumulativeForceAtPT(self,PT):
		f= self.pg.getForceAtPT(PT)
		for s in self.previousPanelGroups:
			f+=s.getForceAtPT(PT)
		return f

	def getCumulativeForceAtRefPT(self):
		return self.getCumulativeForceAtPT(self.pg.refPT)
		
	def setPanelGroup(self,pg):
		self.pg=pg

	def setSection(self,i,s):
		self.section[i]=s

	def setSections(self,s1,s2):
		pass

	def getNode(self,i):
		return self.section[i].orig

	def resetGeometry(self):
		for i in self.section:
			self.section[i].resetGeometry()
		for j in range(0,3):
			self.oVector[j]=self.oVectorInitial[j]

	def getReferenceSection(self):
		if self.integration == '12':
			return self.section[self.section.keys()[-1]]
		else:
			return self.section[0]

	def forces2fem(self):
		ff={}
		f=self.pg.getForceAtRefPT()
		pt=self.pg.refPT
		for m in self.mass:
			f+=m.getForceAtPT(pt)
		if self.integration == '12':
			c=self.ids[1]
		if self.integration=='21':
			c=self.ids[0]
		ff['force']=f
		ff['node']=c

		return ff
		
	def fem2deform(self,U):
		#for i in self.section:
			#self.section[i].moveAll(U[self.ids[i]])
			#self.section[i].defs=copy.copy(U[self.ids[i]])
		if self.integration == '12':
			self.section[1].moveAll(U[self.ids[1]])
			self.section[1].defs=copy.copy(U[self.ids[1]])
			self.section[0].moveAll(U[self.ids[0]])
			self.section[0].defs=copy.copy(U[self.ids[0]])
			self.moveMasses(U[self.ids[1]])
			self.pg.refPT+=U[self.ids[1]][0:3]
			self.orig=self.pg.refPT.copy()

		if self.integration == '21':
			self.section[1].moveAll(U[self.ids[1]])
			self.section[1].defs=copy.copy(U[self.ids[1]])
			
			self.section[0].moveAll(U[self.ids[0]])
			self.section[0].defs=copy.copy(U[self.ids[0]])
			self.moveMasses(U[self.ids[0]])
			self.pg.refPT+=U[self.ids[0]][0:3]
			self.orig=self.pg.refPT.copy()

	def _rotateOVector(self,nodes,U):
		if self.integration == '12':
			i=0
		if self.integration == '21':
			i=1

		rv=U[self.ids[i]][3:6]
		oVector=np.copy(self.oVector)
	
		for j in range(0,3):
			oVector=bc.rotAx(np.array([oVector]),self.csys[j],np.array([0,0,0]),rv[j])[0]

		for j in range(0,3):
			self.oVector[j]=oVector[j]

	def rotateMasses(self,rv):
		csys=np.eye(3)
		for m in self.mass:
			for i in range(0,3):
				tmp=bc.rotAx(np.array([m.p]),csys[i],self.pg.refPT,rv[i])[0]
				for j in range(0,3):
					m.p[j]=tmp[j]

	def moveMasses(self,u):
		self.rotateMasses(u[3:6])
		for m in self.mass:
			m.p+=u[0:3]

	def _rotateOPoint(self,nodes,U):
		if self.integration == '12':
			i=0
		if self.integration == '21':
			i=1

		rv=U[self.ids[i]][3:6]
#		csn=np.copy(self.csysLoc)
		oPoint=np.copy(self.oPoint)
	
		for j in range(0,3):
#			csn=bc.rotAx(csn,self.csys[j],np.zeros(3),rv[j])
			oPoint=bc.rotAx(np.array([oPoint]),self.csys[j],self.orig,rv[j])[0]

		for j in range(0,3):
			self.oPoint[j]=oPoint[j]

	def buildCsys(self):			
		if self.integration == '12':
			i=0
			V0=self.csysPT1-self.orig
		if self.integration == '21':
			i=1
			V0=self.csysPT1-self.orig
		if self.orientation == 'byNode':
			oPoint=self.oPoint
		else:
			oPoint=self.csysPT1+self.oVector
		csn=cs.makeCSYS_from2Vectors(V0, oPoint-self.orig,ax=2)
		for j in range(0,3):
			for k in range(0,3):
				self.csysLoc[j][k]=csn[j][k]
				pass

	