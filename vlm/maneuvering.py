import numpy as np
np.set_printoptions(linewidth=132,precision=9)
import pylab as plt
import copy
import cProfile
from utils.mesh import *
from panelMethod.panelMethod1 import *
from panelMethod.env import GlobalEnvironment
from utils.b3Vect import B3vect 
from utils.bCad import BCad
from resampleFoil import resample
from utils.bCad import *
from am_analysis.amFunctions import *
from utils.panukl2pm import *

from fem.fem4 import FEM1
from fem.beam3d_nl import Beam3D
from utils.interactor import PanelBeamInteractor
from utils.loftedSurface import *
from utils.elasticWing import ElasticWing1
import utils.csys as cs
from utils.beamCalc import BeamCrossSection as BeamCalc
from fem.utils.writeBDF import BDFWriter
from utils.simpleSTEPReader import extractPointsFromSTEPFile
from utils.geometryGen import *
from flightLoads.pointLoads import *
from utils.stl_io import extractPtsFromSTL


import affinity
import multiprocessing
 
affinity.set_process_affinity_mask(0,2**multiprocessing.cpu_count()-1)


from utils.bCad import *
from inertia.inertiaModel1 import InertiaModel
from utils.transform import rotateVectorByAngles
from panelMethod.doubletSourcePanels import doubletSourcePanel
from panelMethod.graphics import *
from panelMethod.panelMethod1 import *
from utils.panukl2pm import *

from am_analysis.geometry.am008.am008_makeGeometry_maneuver1 import makeGeometry

class FlightTimeStep(object):
	def __init__(self):
		self.geometryFile=''
		self.t=0
		self.dt=0.0
		self.CG=np.zeros(3)
		self.COR=np.zeros(3)
		self.G=np.array([0,0,-1])
		self.THRUST=np.array([1,0,0])
		self.POSITION=np.zeros(6)
		self.DISPLACEMENT=np.zeros(6)
		self.FORCE=np.zeros(6)
		self.FLIGHTVEL=np.zeros(3)
		self.dVEL=np.zeros(6)
		self.ACC=np.zeros(6)
		self.V=np.zeros(3)
		self.OMEGA=np.zeros(3)
		self.panels={}
		self.pg=PanelGroup()
		
	def setAircraft(self,pr1,im):#,aircraft,im):
		self.pr1=pr1
		self.aircraft=pr1.panelGroups['AIRCRAFT']
		self.im=im
	
	def compute(self,nextStep=None):

		self.pr1.setFreeVelocity(self.V)
		self.pr1.dom1.COR=self.COR
		#self.pr1.dom1.omega=self.OMEGA
		self.pr1.getSolution()

		self.FORCE=self.aircraft.getForceAtPT(self.CG)
		self.FORCE[0:3]+=self.G
		self.FORCE[0:3]+=self.THRUST
		
		self.ACC[0:3]=self.FORCE[0:3]/960.#self.im.getMass()
		self.ACC[3:6]=self.FORCE[3:6]/5.3e2
		
		self.dVEL=self.ACC*self.dt
		self.FLIGHTVEL+=self.dVEL[0:3]
		self.DISPLACEMENT=self.ACC*0.5*self.dt**2
		self.POSITION+=self.DISPLACEMENT
		self.POSITION[0:3]+=self.FLIGHTVEL*self.dt
		
		cargs=['gid','vloc','forceVector','pressureCoef','pressure','S','G']
		for i in self.aircraft.panels:
			p=self.aircraft.panels[i]
			if not p.isWake:
				pp=doubletSourcePanel(p.pts[0],p.pts[1],p.pts[2],p.pts[3])
				pp.pts=rotateX(pp.pts,self.CG,self.POSITION[3])
				pp.pts=rotateY(pp.pts,self.CG,self.POSITION[4])
				pp.pts=rotateZ(pp.pts,self.CG,self.POSITION[5])
				
				pp.pts+=self.POSITION[0:3]
				for ar in cargs:
					setattr(pp,ar,getattr(p,ar))
				self.pg.addPanel(pp)
		
		if nextStep is not None:
			nextStep.V=np.copy(rotateVectorByAngles(self.V,self.DISPLACEMENT[3:6]))
			nextStep.FLIGHTVEL=np.copy(rotateVectorByAngles(self.FLIGHTVEL,-self.DISPLACEMENT[3:6]))
			nextStep.G=np.copy(rotateVectorByAngles(self.G,-self.DISPLACEMENT[3:6]))
			#nextStep.OMEGA=self.OMEGA+self.VEL[3:6])
			#nextStep.OMEGA=np.copy(self.VEL[3:6])
			nextStep.POSITION=np.copy(self.POSITION)
	
	def loadAero(self,fname):
		self.pr1.loadAeroFromFile(fname)

class Maneuver1(object):
	def __init__(self):
		self.basename=''
		self.time={}
		self.control={}
		self.CG=np.zeros(3)
		self.COR=np.zeros(3)
		self.panels={}
		self.im=InertiaModel()
		
	def setAeroDomain(self,pr1):	
		self.pr1=pr1
		
	def initSteps(self,v0,fVel,td,CG,THRUST):
		self.control=td
		for i in td:
			self.time[i]=FlightTimeStep()
			self.time[i].V=np.copy(v0)
			self.time[i].FLIGHTVEL=np.copy(fVel)
			self.time[i].CG=np.copy(CG)
			self.time[i].THRUST=np.copy(THRUST)
			self.time[i].dt=td[i]['dt']
			self.time[i].geometryFile=td[i]['geometryFile']

	
	def compute(self,pr1,plot=False):
		pr1.__init__()
		ts=self.time[0]
		
		pr1.setFreeVelocity(ts.V)
		
		
		makeGeometry(pr1,[ts.geometryFile])
		ts.setAircraft(pr1,None)#,pr1.AIRCRAFT,None)
		self.time[0].compute(nextStep=self.time[1])
		lastGeometry=ts.geometryFile
		pr1.saveAeroToFile(self.basename+str(0)+'.aero')

		if plot:
			i=0
			PGtoVTK(self.basename+str(i)+'.vtm',self.time[i].pr1.panelGroups,plotWake=False)
		del(self.time[i].pg)
		self.time[i].pg=None
		for i in range(1,len(self.time)-1):
			ts=self.time[i]
			if lastGeometry != ts.geometryFile:
				pr1.__init__()
				
				makeGeometry(pr1,[ts.geometryFile])
			pr1.setFreeVelocity(ts.V)
			ts.setAircraft(pr1,None)#,pr1.AIRCRAFT,None)
			self.time[i].compute(nextStep=self.time[i+1])
			pr1.saveAeroToFile(self.basename+str(i)+'.aero')
			if plot:
				PGtoVTK(self.basename+str(i)+'.vtm',{'AIRCRAFT':self.time[i].pg},plotWake=False)
			del(self.time[i].pg)
			self.time[i].pg=None			
			
	def plot(self,basename=''):
		for i in self.time:
			PGtoVTK(basename+str(i)+'.vtm',self.time[i].pr1.panelGroups,plotWake=False)

	def plotPost(self):
		pr1.__init__()
		ts=self.time[0]
		
		pr1.setFreeVelocity(ts.V)
		
		for i in range(0,len(self.time)-1):
			ts=self.time[i]
			pr1.__init__()
			makeGeometry(pr1,[ts.geometryFile])
			pr1.setFreeVelocity(ts.V)
			ts.setAircraft(pr1,None)#,pr1.AIRCRAFT,None)
			pr1.loadAeroFromFile(self.basename+str(i)+'.aero')
			PGtoVTK(self.basename+str(i)+'.vtm',{'AIRCRAFT':self.time[i].pg},plotWake=False)
		
pr1=VlmProblem()
m=Maneuver1()
td={}

for i in range(0,30):
	td[i]={}
	td[i]['dt']=0.04
	td[i]['geometryFile']='am008_incidence000_elevator030up.dat'

for i in range(30,1000):
	td[i]={}
	td[i]['dt']=0.04
	td[i]['geometryFile']='am008_incidence000_elevator030up.dat'


td[0]['geometryFile']='am008_incidence000_elevator000up.dat'
td[1]['geometryFile']='am008_incidence000_elevator005up.dat'
td[2]['geometryFile']='am008_incidence000_elevator010up.dat'
td[3]['geometryFile']='am008_incidence000_elevator015up.dat'
td[4]['geometryFile']='am008_incidence000_elevator020up.dat'
td[5]['geometryFile']='am008_incidence000_elevator025up.dat'
td[6]['geometryFile']='am008_incidence000_elevator030up.dat'


td[30]['geometryFile']='am008_incidence000_elevator030up.dat'
td[31]['geometryFile']='am008_incidence000_elevator025up.dat'
td[32]['geometryFile']='am008_incidence000_elevator020up.dat'
td[33]['geometryFile']='am008_incidence000_elevator015up.dat'
td[34]['geometryFile']='am008_incidence000_elevator010up.dat'
td[35]['geometryFile']='am008_incidence000_elevator005up.dat'
td[36]['geometryFile']='am008_incidence000_elevator000up.dat'






m.basename='am_analysis/data/am008_loop/x'
m.initSteps(np.array([100.0,0,0]),np.array([-100.0,0,0]),td,np.array([2.8,0,1.06]),np.array([-2700.0,0,0]))
#m.compute(pr1,plot=False)
m.compute(pr1,plot=True)
#m.plot('xxxxxx_')
m.plotPost()