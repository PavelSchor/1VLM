import sys
#sys.path.insert(0, '/usr/local/lib/python2.7/site-packages/')
import numpy as np
import utils.csys as cs
from collections import OrderedDict
import json_tricks.np as json# import dump, dumps, load, loads, strip_comments
from utils.transform import rotateVectorByAngles
	
class FlightModel1(object):
	def __init__(self):
		self.postProcess=True
		self.solutionDict={'solAero':'linear'}
		self.mode='elastic'
		self.geometryFile=''
		self.t=0
		self.dt=0.0
		self.CSYS_BODY=np.eye(3)
		self.CSYS_PANELS=np.array([ [-1,0,0],[0,-1,0],[0,0,1] ])
		self.CSYS_WIND=np.array([ [-1,0,0],[0,-1,0],[0,0,1] ])

		self.CSYS_EARTH=np.eye(3)
		self.CG=np.zeros(3)
		self.COR=np.zeros(3)
		self.G=np.array([0,0,-1])
		self.THRUST=np.array([-1,0,0])*0
		self.THRUST_PT=np.array([6.587,0.,1.501])

		self.POSITION=np.zeros(6)
		self.DISPLACEMENT=np.zeros(6)
		self.FORCE=np.zeros(6)
		self.FLIGHTVEL=np.zeros(3)
		self.dVEL=np.zeros(6)
		self.ACC=np.zeros(6)
		self.V=np.zeros(3)
		
		self.MASS=500.
		self.panels={}
		self.ACC_INERT=np.zeros(3)
		self.ACC_EARTH=np.zeros(3)
		self.VEL_EARTH=np.zeros(3)
		self.POS_EARTH=np.zeros(3)
		self.OMEGA=np.zeros(3)
		self.RPOS_EARTH=np.zeros(3)

		self.wakeDict={}
		
		self.elevator=0.0
		self.II=np.eye(3)*3000
		self.camPosition=np.array([0,0,0])
		self.camFocalPoint=np.array([50,0,0])
		self.camUpview=np.array([0,0,1])
		
		self.camGlobalPosition=np.zeros(3)
		self.camGlobalFocalPoint=np.zeros(3)
		self.camGlobalUpview=np.zeros(3)
		self.restrictedDOF=np.zeros(6,dtype=np.bool)
		
		self.joystick=None
		
	def getGravity(self,csys='body'):
		if csys =='body':
			return cs.transformVector(np.array([0,0,-9.80665]),self.CSYS_BODY,self.CSYS_EARTH)
		elif csys =='panels':
			return cs.transformVector(np.array([0,0,-9.80665]),self.CSYS_PANELS,self.CSYS_EARTH)
			
	def assembleStateVector(self):
			
		#		x(1)    = 		Body-axis x inertial velocity, ub, m/s
		#		x(2)    =		Body-axis y inertial velocity, vb, m/s
		#		x(3)    =		Body-axis z inertial velocity, wb, m/s
		#		x(4)    =		North position of center of mass WRT Earth, xe, m
		#		x(5)    =		Eaself position of center of mass WRT Earth, ye, m
		#		x(6)    =		Negative of c.m. altitude WRT Earth, ze = -h, m
		#		x(7)    =		Body-axis roll rate, pr, rad/s
		#		x(8)    =		Body-axis pitch rate, qr, rad/s
		#		x(9)    =		Body-axis yaw rate, rr,rad/s
		#		x(10)   =		Roll angle of body WRT Earth, phir, rad
		#		x(11)   =		Pitch angle of body WRT Earth, thetar, rad
		#		x(12)   =		Yaw angle of body WRT Earth, psir, rad
		self.X=np.zeros(12)
		
		self.X[0:3]=np.copy(self.V)
		self.X[3:6]=np.copy(self.POS_EARTH)
		self.X[6:9]=np.copy(self.OMEGA)
		self.X[9:12]=np.copy(self.RPOS_EARTH)
		
		
		#self.FORCE_AERO=self.aircraft.getForceAtPT(self.CG)
		#vd=self.pr1.getFreeVelocity()
		#self.FORCE_AERO_VISC=(0.5*self.pr1.getRho()*np.linalg.norm(vd)**2)*(vd/np.linalg.norm(vd))*17.44*0.0011
		#self.FORCE=self.aircraft.getForceAtPT(self.CG)
		self.FORCE[0:3]=np.copy(self.dom1.getForce())
		self.FORCE[3:6]=np.copy(self.dom1.getMoment())
		self.FORCE_AERO=self.FORCE.copy()

		self.FORCE[0:3]=cs.transformVector(self.FORCE[0:3],self.CSYS_BODY,self.CSYS_PANELS)
		#self.FORCE[0:3]+=cs.transformVector(self.FORCE_AERO_VISC,self.CSYS_BODY,self.CSYS_PANELS)
		
		self.FORCE[3:6]=cs.transformVector(self.FORCE[3:6],self.CSYS_BODY,self.CSYS_PANELS)
		self.G=self.getGravity(csys='body')*self.MASS
		self.FORCE[0:3]+=self.G
		self.FORCE[0:3]+=cs.transformVector(self.THRUST,self.CSYS_BODY,self.CSYS_PANELS)
		
		rThr=self.THRUST_PT-self.CG
		
		self.FORCE[3:6]+=cs.transformVector(np.cross(rThr,self.THRUST),self.CSYS_BODY,self.CSYS_PANELS)
		
		#self.FORCE[self.restrictedDOF]=0
		
		

		def skew(x):#antisymmetric matrix
			return np.array([[  0,-x[2], x[1]],
					[ x[2], 0, -x[0]],
					[-x[1], x[0], 0]])
		self.ACC_LIN=self.FORCE[0:3]/self.MASS - np.dot(skew(self.OMEGA),self.V)

		self.ACC_ROT=np.dot(np.linalg.inv(self.II), self.FORCE[3:6] - np.dot( np.dot(skew(self.OMEGA),self.II),self.OMEGA))
		#self.ACC_ROT=self.FORCE[3:6]/5.3e3
		#self.ACC_ROT[[0,2]]=0.0

		self.ACC_LIN[self.restrictedDOF[0:3]]=0
		self.ACC_ROT[self.restrictedDOF[3:6]]=0
			   
		self.XD=np.zeros(12)
		self.XD[0:3]=self.ACC_LIN
		self.XD[3:6]=self.V
		self.XD[6:9]=self.ACC_ROT
		self.XD[9:12]=self.OMEGA


		return self.XD

	def deflectControls(self):
		#self.controls.pr1=self.pr1
		#l=['Aileron','Elevator','Rudder']
		#for i in l:
			#if not getattr(self.controls,'limiter'+i) is None:
				#lim=getattr(self.controls,'limiter'+i)
				#lim.dt=self.dt
				#lim.t=self.t
		self.controls.setDeflectionsInTime(self.t)

	def compute(self):
		
		self.dom1.setFreeVelocity(cs.transformVector(self.V,self.CSYS_PANELS,self.CSYS_BODY)*-1)
		self.dom1.COR=np.copy(self.CG)
		self.dom1.OMEGA= cs.transformVector(self.OMEGA,self.CSYS_PANELS,self.CSYS_BODY)*-1.

		if self.joystick is not None:
			self.joystick.applyInputs(self.dom1.viscCL_correction)
			
		self.dom1.compute()

		self.assembleStateVector()
	
		DX=self.dt*self.XD
		XNEW=self.X+DX
		self.V_OLD=np.copy(self.X[0:3])
		self.V[0]=XNEW[0];self.V[1]=XNEW[1];self.V[2]=XNEW[2]
		self.OMEGA=np.copy(XNEW[6:9])
		
		self.POS_EARTH=np.copy(  self.POS_EARTH+cs.transformVector(DX[3:6],self.CSYS_EARTH,self.CSYS_BODY) )
		

		
		
		self.RPOS_EARTH= np.copy(  self.RPOS_EARTH+cs.transformVector(DX[9:12],self.CSYS_EARTH,self.CSYS_BODY))

		#self.CSYS_EARTH=cs.rotateCSYSByAngles(self.CSYS_EARTH,cs.transformVector(DX[9:12],self.CSYS_BODY,self.CSYS_EARTH))
		self.CSYS_EARTH=np.copy( cs.rotateCSYSByAngles(self.CSYS_EARTH,DX[9:12]*-1)) 
		self.CSYS_EARTH=cs.grammSmithOrthoNormalize(self.CSYS_EARTH)
		self.CSYS_WIND=np.copy(cs.rotateCSYSByAngles(self.CSYS_WIND,DX[9:12]))
		self.CSYS_WIND=cs.grammSmithOrthoNormalize(self.CSYS_WIND)



	
		self.ACC_EARTH=np.copy( (cs.transformVector(self.V,self.CSYS_EARTH,self.CSYS_BODY)- cs.transformVector(self.V_OLD,self.CSYS_EARTH,self.CSYS_BODY))/self.dt)
		self.ACC_INERT=np.copy( cs.transformVector(self.ACC_EARTH,self.CSYS_BODY,self.CSYS_EARTH) )

		self.ACC_INERT=np.copy( cs.transformVector(self.ACC_EARTH,self.CSYS_BODY,self.CSYS_EARTH) )


		self.camGlobalPosition=cs.transformVector(self.camPosition,self.CSYS_EARTH,self.CSYS_BODY)+self.POS_EARTH
		self.camGlobalFocalPoint=cs.transformVector(self.camFocalPoint,self.CSYS_EARTH,self.CSYS_BODY)+self.POS_EARTH
		self.camGlobalUpview=cs.transformVector(self.camUpview,self.CSYS_EARTH,self.CSYS_BODY)
		

			

	def getCameraPos(self):
		return self.camGlobalPosition
	
	def getCameraFocal(self):
		return self.camGlobalFocalPoint
	
	def getCameraViewUp(self):
		return self.camGlobalUpview

		
		
	def saveJSON(self,fname):
		self.gravity_s=self.getGravity(csys='panels')
		self.acceleration_s=cs.transformVector(self.ACC_EARTH,self.CSYS_PANELS,self.CSYS_EARTH)
		self.angularAcceleration_s=cs.transformVector(self.ACC_ROT,self.CSYS_PANELS,self.CSYS_BODY)
		self.angularVelocity_s=cs.transformVector(self.OMEGA,self.CSYS_PANELS,self.CSYS_BODY)
		self.controlInputs=self.controls.getDeflectionDict()
		vvv=self.dom1.getFreeVelocity()
		self.alpha=np.arctan2(vvv[2],np.abs(vvv[0]))
		self.beta=np.arctan2(vvv[1],np.abs(vvv[0]))
		self.velocity=np.linalg.norm(vvv)
		req=['t','velocity','alpha','beta','geometryFile','controlInputs',"CSYS_AERO","CSYS_BODY","CSYS_EARTH","CG","V","OMEGA","FORCE","G","POS_EARTH","RPOS_EARTH","X","XD",'ACC_INERT','ACC_EARTH','gravity_s','acceleration_s','angularAcceleration_s','angularVelocity_s','FORCE_AERO','THRUST','ACC_ROT']
		reqd=OrderedDict()
		for i in req:
			try:
				reqd[i]=getattr(self,i)
			except:
				pass
		for i in self.dom1.panelGroups:
			reqd[f'force_{i}']=self.dom1.panelGroups[i].getForceAtRefPT()

		
		with open(fname,'w') as data_file:
			data_file.write(json.dumps(reqd,indent=4))
			
	def loadJSON(self,fname):
		with open(fname,'r') as f:
			d=json.load(f)
		for i in d:
			setattr(self,i,d[i])


