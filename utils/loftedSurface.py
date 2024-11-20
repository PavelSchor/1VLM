import numpy as np
from utils.b3Vect import B3vect
import utils.bCad
from utils.bCad import *
import copy
bc=utils.bCad.BCad()
vect=B3vect()


class LoftSection(object):
#For aeroelastic purposes
	def __init__(self):
		self.allowDeform=True
		self.points={}
		self.pointsSurf={}
		self.initialPoints={}
		self.initialPointsSurf={}
		self.orig=np.zeros(3)
		self.oPT=np.zeros(3)
		self.npts=0
		self.nptsSurf=0
		self.csys=np.array([ [1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
		pass

	def addPoint(self,pt):
		if pt.shape==(3,):
			self.points[self.npts]=pt
			self.initialPoints[self.npts]=copy.copy(pt)
			self.npts+=1
		else:
			for i in range(0,pt.shape[0]):
				self.addPoint(pt[i])

	def addSurfacePoint(self,pt):
		if pt.shape==(3,):
			self.pointsSurf[self.nptsSurf]=pt
			self.initialPointsSurf[self.nptsSurf]=copy.copy(pt)
			self.nptsSurf+=1
		else:
			for i in range(0,pt.shape[0]):
				self.addSurfacePoint(pt[i])


	def pts2arr(self):
		arr=np.zeros((self.npts,3))
		for i in range(0,self.npts):
			arr[i]=self.points[i]
		return arr


	def ptsSurf2arr(self):
		arr=np.zeros((self.nptsSurf,3))
		for i in range(0,self.nptsSurf):
			arr[i]=self.pointsSurf[i]
		return arr


	def arr2pts(self,arr):
		if arr.shape[0]==self.npts:
			for i in range(0,self.npts):
				self.points[i]=arr[i]
		else:
			self.npts=0
			self.addPoint(arr)

	def rotate(self,rv):
		csys=np.eye(3)
		for i in range(0,3):
			
			tmp=bc.rotAx(self.pts2arr(),csys[i],self.orig,rv[i])
			self.tmp=tmp
		
			for j in range(0,self.npts):
				for k in range(0,3):
					self.points[j][k]=tmp[j][k]

	def translate(self,tv):
		for j in range(0,self.npts):
				self.points[j]+=tv
		#self.points+=tv
		
	def move(self,v):
		self.rotate(v[3:6])
		self.translate(v[0:3])
		

	def moveAll(self,v):
		#self.move(v)
		self.rotate(v[3:6])
		self.orig+=v[0:3]
		self.translate(v[0:3])
		
	def getArea(self):
		#return bc.polyArea(self.ptsSurf2arr())
		return polyArea(self.ptsSurf2arr()[:,[0,2]])
	
	def getCircumference(self):
		return bc.polyCircumference(self.ptsSurf2arr())

	def getJk(self,t):
		U=self.getArea()
		s=self.getCircumference()
		return (4.*U**2.0)/(s/t)
	
	def resetGeometry(self):
		for i in range(0,self.npts):
			for j in range(0,3):
				self.points[i][j]=self.initialPoints[i][j]


class TwoSectionsLoft(object):

	def __init__(self):
		self.nInterpolants=0
		pass

	def setSectionA(self,A):
		self.sectionA=A

	def setSectionB(self,B):
		self.sectionB=B

	def setSections(self,A,B):
		self.setSectionA(A)
		self.setSectionB(B)

	def setNumberOfIterpolants(self,n):
		self.nInterpolants=n

	def getInterpolants(self):
		res={}
		ni=2+self.nInterpolants

		for i in range(1,ni-1):
			res[i]=np.zeros(self.sectionA.shape)
			for j in range(0,self.sectionA.shape[0]):
				v=self.sectionB[j]-self.sectionA[j]
				res[i][j]=self.sectionA[j]+v*(float(i)/float(ni-1))
#				print 'i, ni, i.ni',i, ni, float(i)/float(ni-1)
		res[0]=self.sectionA
		res[ni-1]=self.sectionB
		return res

class MultiSectionSurface(object):

	def __init__(self):
		self.sections={}
		self.nSec=0
		self.rudder=None
		pass


	def addSection(self,A,B,n=0):
		self.sections[self.nSec]=TwoSectionsLoft()
		self.sections[self.nSec].setSections(A,B)
		self.sections[self.nSec].setNumberOfIterpolants(n)
		self.nSec+=1

	def getRawInterpolants(self):
		res={}
		ni=0
		sres=self.sections[0].getInterpolants()
		for j in range(0, len(sres)):
			res[ni]=sres[j]
			ni+=1
		for i in range(1,self.nSec):
			sres=self.sections[i].getInterpolants()
			for j in range(0, len(sres)):
				resTmp=sres[j]
				diff=np.abs(resTmp-res[ni-1])
#				print 'diff.max() ', diff.max()
				if diff.max() > 1e-10:
					res[ni]=resTmp
					ni+=1
#		res[ni]=sres[j+1]
#		ni+=1

		return res

	def getInterpolants(self):
		if self.rudder is None:
			return self.getRawInterpolants()
		else:
			r={}
			res= self.getRawInterpolants()
			self.rudder.deflect(res,False)
			for i in range(0,len(res)):
				r[i]=res[i][self.rudder.activePts]
			
			return r