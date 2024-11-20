import numpy as np
from pyNastran.bdf import bdf
from scipy.spatial import cKDTree
import collections
import vtk
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.field_writer_double import print_card_double

import affinity
import multiprocessing
affinity.set_process_affinity_mask(0,2**multiprocessing.cpu_count()-1)


def uniqueRows(a):
	a = np.ascontiguousarray(a)
	unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
	return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))

class CSHELL(object):
	def __init__(self,ids):
		self.ids=ids
		pass

	def getCPoint(self,grid):
		pts=uniqueRows(grid[self.ids])
		return np.sum(pts,axis=0)/len(pts)


class BDFSHELLSHANDLER(object):
	def __init__(self):
		pass


	def loadInput(self,fname):
		
		self.bdfIn=bdf.BDF()
		self.bdfIn.read_bdf(fname)
		
		
		self.nodeHT={}
		self.nodeHTinv={}
		j=0
		for i in self.bdfIn.nodes:
			self.nodeHT[j]=self.bdfIn.nodes[i].nid
			self.nodeHTinv[self.bdfIn.nodes[i].nid]=j
			j+=1
		
		self.elHT={}
		self.elHTinv={}
		elK=self.bdfIn.elements.keys()
		j=0
		for i in range(0,len(elK)):
			e=self.bdfIn.elements[elK[i]]
			if e.type=='CQUAD4' or e.type=='CTRIA3':
				self.elHT[j]=e.eid
				self.elHTinv[e.eid]=j
				j+=1
		
		for i in self.elHT.values():
			e=self.bdfIn.elements[i]
			e.pressureLoad=0.0
			e.pressureCoef=1.0	
		

	def mapPLOAD(self,pr1,scaleGeo=1.,scalePress=1.,k=10,pp=2,set0=None,mn=(2,2),method='centroids'):
		for i in self.elHT.values():
			e=self.bdfIn.elements[i]
			e.pressureLoad=0.0
			e.pressureCoef=0.0	


		#centroidsIn=np.zeros((len(self.elHT),3))
		#for i in range(0,len(self.elHT)):
			#centroidsIn[i]=self.bdfIn.elements[self.elHT[i]].Centroid()

		if type(set0) is list:
			self.mappedEl=list(set(set0).intersection( set(self.elHT.values() )))
		elif type(set0) is int:
			self.mappedEl=list(set(self.bdfIn.sets[set0].ids).intersection(set(self.elHT.values() )))
		else:
			self.mappedEl=self.elHT.values()

		
		p1=pr1.getBodyPanels()
		nmp=mn[0]*mn[1]
		nco=len(p1)*nmp
		g1=np.zeros((nco,3))
		l1=np.zeros(nco)
		c1=np.zeros(nco)

		for i in p1:
			g1[i*nmp:(i+1)*nmp]=p1[i].getSubCentroids(mn[0],mn[1])
			l1[i*nmp:(i+1)*nmp]=np.ones(nmp)*p1[i].getPressureLoad()
			c1[i*nmp:(i+1)*nmp]=np.ones(nmp)*p1[i].getPressureCoef()
			#g1[i]=p1[i].getCpoint()
			#l1[i]=p1[i].getPressureLoad()
			#c1[i]=p1[i].getPressureCoef()
			
		self.mappedPoints=np.zeros((len(self.mappedEl),3))
		self.mappedCP=np.zeros(len(self.mappedEl))
		self.mappedPLOAD=np.zeros(len(self.mappedEl))

		if method=='centroids':
			for i in range(0,len(self.mappedEl)):
				e=self.bdfIn.elements[self.mappedEl[i]]
				self.mappedPoints[i]=e.Centroid()


			self.mappedPoints*=scaleGeo

			tree = cKDTree(g1)


			d, inds = tree.query(self.mappedPoints, k =k,p=pp)
			w = 1.0 / d**2
			if k==1:
				self.mappedPLOAD = l1[inds]
				self.mappedCP = c1[inds]
			else:
				self.mappedPLOAD = np.sum(w * l1[inds], axis=1) / np.sum(w, axis=1)
				self.mappedCP = np.sum(w * c1[inds], axis=1) / np.sum(w, axis=1)
			self.mappedPLOAD*=scalePress
			
			for i in range(0,len(self.mappedEl)):
				e=self.bdfIn.elements[self.mappedEl[i]]
				e.pressureCoef=self.mappedCP[i]
				e.pressureLoad=self.mappedPLOAD[i]
			
		if method=='nodes':
			isCTRIA3 = np.zeros(len(self.mappedEl),dtype=np.bool)
			isCQUAD4 = np.zeros(len(self.mappedEl),dtype=np.bool)
			for i in range(0,len(self.mappedEl)):
				e=self.bdfIn.elements[self.mappedEl[i]]
				if e.type=='CQUAD4':
					isCQUAD4[i]=True
				if e.type=='CTRIA3':
					isCTRIA3[i]=True
			
			self.mappedPoints=np.zeros((isCQUAD4.sum()*4+isCTRIA3.sum()*3,3))
			iHead=0
			for i in range(0,len(self.mappedEl)):
				e=self.bdfIn.elements[self.mappedEl[i]]
				for j in range(0,len(e.nodes)):
					self.mappedPoints[iHead]=e.nodes[j].xyz
					iHead+=1

			self.mappedPoints*=scaleGeo
			tree = cKDTree(g1)
			d, inds = tree.query(self.mappedPoints, k =k,p=pp)
			w = 1.0 / d**2
			if k==1:
				self.mappedPLOAD = l1[inds]
				self.mappedCP = c1[inds]
			else:
				self.mappedPLOAD = np.sum(w * l1[inds], axis=1) / np.sum(w, axis=1)
				self.mappedCP = np.sum(w * c1[inds], axis=1) / np.sum(w, axis=1)
				
			self.mappedPLOAD*=scalePress
			iHead=0
			for i in range(0,len(self.mappedEl)):
				e=self.bdfIn.elements[self.mappedEl[i]]
				if isCTRIA3[i]:
					e.pressureCoef=self.mappedCP[iHead:iHead+3].sum()/3.
					e.pressureLoad=self.mappedPLOAD[iHead:iHead+3].sum()/3.
					iHead+=3
				if isCQUAD4[i]:
					e.pressureCoef=self.mappedCP[iHead:iHead+4].sum()/4.
					e.pressureLoad=self.mappedPLOAD[iHead:iHead+4].sum()/4.
					iHead+=4

	def pload2ToFile(self,fname,sid=2):
		f=open(fname,'w')
		f.write(print_card_8(['LOAD', sid, 1., 1., 1]))
		for i in self.mappedEl:
			e=self.bdfIn.elements[i]
			f.write(print_card_8(['PLOAD2',sid ,e.pressureLoad,e.eid]))
		f.close()
		
	def toVTK(self,fname):
		def getNodePos(ii):
			return self.nodeHTinv[ii]
		
		
		points = vtk.vtkPoints()
		cells = vtk.vtkCellArray()

		for i in self.bdfIn.nodes:
			ids=points.InsertNextPoint(self.bdfIn.nodes[i].xyz  )

		arr_cp=vtk.vtkFloatArray()
		arr_cp.SetName('CP')
		arr_pload=vtk.vtkFloatArray()
		arr_pload.SetName('PLOAD')
		arr_eids=vtk.vtkFloatArray()
		arr_eids.SetName('EID')

		#for i in self.elHT:
			#e=self.bdfIn.elements[self.elHT[i]]
		for i in self.mappedEl:
			e=self.bdfIn.elements[i]
			if e.type=='CQUAD4':
				quad=vtk.vtkQuad()
				quad.GetPointIds().SetId(0,getNodePos(e.node_ids[0]))
				quad.GetPointIds().SetId(1,getNodePos(e.node_ids[1]))
				quad.GetPointIds().SetId(2,getNodePos(e.node_ids[2]))
				quad.GetPointIds().SetId(3,getNodePos(e.node_ids[3]))
				cells.InsertNextCell(quad)
			if e.type=='CTRIA3':
				triangle = vtk.vtkTriangle();

				triangle.GetPointIds().SetId(0,getNodePos(e.node_ids[0]))
				triangle.GetPointIds().SetId(1,getNodePos(e.node_ids[1]))
				triangle.GetPointIds().SetId(2,getNodePos(e.node_ids[2]))
				cells.InsertNextCell(triangle);

			
			arr_eids.InsertNextValue(e.eid)
			arr_cp.InsertNextValue(e.pressureCoef)
			arr_pload.InsertNextValue(e.pressureLoad)
		

			
		polydata = vtk.vtkPolyData()
		polydata.SetPoints(points)
		polydata.SetPolys(cells)
		#polydata.Modified()

		
		polydata.GetCellData().AddArray(arr_eids)
		polydata.GetCellData().AddArray(arr_cp)
		polydata.GetCellData().AddArray(arr_pload)
		
		#polydata.Modified()
		if vtk.VTK_MAJOR_VERSION <= 5:
			polydata.Update()
		
		writer = vtk.vtkXMLPolyDataWriter();
		writer.SetFileName(fname);
		if vtk.VTK_MAJOR_VERSION <= 5:
			writer.SetInput(polydata)
		else:
			writer.SetInputData(polydata)
		writer.Write()







		
		
		
		