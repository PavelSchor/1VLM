from utils.mesh import *
from vlm.vlmMethod1 import *
from flight.fcs1 import FlightControlSystem1
from flight.fdm1 import FlightModel1
import numpy as np
np.setbufsize(8192*2)
import pylab as plt
import time
from tqdm import tqdm
from scipy.optimize import root




cg=np.array([0.25,0,0])
chord=1.
span=10.
alpha=np.deg2rad(10)
beta=np.deg2rad(0)
sp=35.
q=0.5*1*sp**2
VEL=np.array([np.cos(alpha),-np.sin(beta),np.sin(alpha)])*sp
pr1=VlmProblem()
self=pr1.dom1
pr1.dom1.COR=cg

xDir=np.array([1,0,0])
yDir=np.array([0,1,0])
zDir=np.array([0,0,1])

wing_l=PanelGroup()
nc=21
xf=np.linspace(0,1,nc)
pts1=np.zeros((nc,3)); pts1[:,0]=np.linspace(0,1,nc);pts1[:,2]=0.05*np.sin(xf*np.pi); pts1[:,1]=-5.
pts2=np.zeros((nc,3)); pts2[:,0]=np.linspace(0,1,nc);pts2[:,2]=0.05*np.sin(xf*np.pi);pts2[:,1]=-0.5
ll=[pts1,pts2]
nip=[10]
makeMeshFromArrayList(pr1,wing_l,ll,nip=nip)
appendLoftedWakeRegKnownLasts(pr1,wing_l,1,wakeVector=np.array([30,0,0]))

wing_r=PanelGroup()
nc=21
pts1=np.zeros((nc,3)); pts1[:,0]=np.linspace(0,1,nc);pts1[:,2]=0.05*np.sin(xf*np.pi);pts1[:,1]=0.5
pts2=np.zeros((nc,3)); pts2[:,0]=np.linspace(0,1,nc);pts2[:,2]=0.05*np.sin(xf*np.pi);pts2[:,1]=5.
ll=[pts1,pts2]
nip=[10]
makeMeshFromArrayList(pr1,wing_r,ll,nip=nip)
appendLoftedWakeRegKnownLasts(pr1,wing_r,1,wakeVector=np.array([30,0,0]))

htu_l=PanelGroup()
nc=11
pts1=np.zeros((nc,3)); pts1[:,0]=7+np.linspace(0,0.5,nc);pts1[:,2]=2.5;pts1[:,1]=-2
pts2=np.zeros((nc,3)); pts2[:,0]=7+np.linspace(0,0.5,nc);pts2[:,2]=2.5;pts2[:,1]=0.
ll=[pts1,pts2]
nip=[10]
makeMeshFromArrayList(pr1,htu_l,ll,nip=nip)
appendLoftedWakeRegKnownLasts(pr1,htu_l,1,wakeVector=np.array([30,0,0]))

htu_r=PanelGroup()
nc=11
pts1=np.zeros((nc,3)); pts1[:,0]=7+np.linspace(0,0.5,nc);pts1[:,2]=2.5;pts1[:,1]=0.
pts2=np.zeros((nc,3)); pts2[:,0]=7+np.linspace(0,0.5,nc);pts2[:,2]=2.5;pts2[:,1]=2
ll=[pts1,pts2]
nip=[10]
makeMeshFromArrayList(pr1,htu_r,ll,nip=nip)
appendLoftedWakeRegKnownLasts(pr1,htu_r,1,wakeVector=np.array([30,0,0]))


vtu=PanelGroup()
nc=11
pts1=np.zeros((nc,3)); pts1[:,0]=7+np.linspace(0,0.5,nc);pts1[:,2]=0.5;
pts2=np.zeros((nc,3)); pts2[:,0]=7+np.linspace(0,0.5,nc);pts2[:,2]=2.5;
ll=[pts2,pts1]
nip=[10]
makeMeshFromArrayList(pr1,vtu,ll,nip=nip)
appendLoftedWakeRegKnownLasts(pr1,vtu,1,wakeVector=np.array([30,0,0]))

nc=221
xf=np.linspace(0,1,nc)
fuselage_v=PanelGroup()

pts1=np.zeros((nc,3)); pts1[:,0]=xf*11-3.5;pts1[:,2]=0.5*np.sin(xf*np.pi)+0.2; pts1[:,2][pts1[:,0]>=7]=0.5;
pts2=np.zeros((nc,3)); pts2[:,0]=xf*11-3.5;pts2[:,2]=0.
pts3=np.zeros((nc,3)); pts3[:,0]=xf*11-3.5;pts3[:,2]=-0.5*np.sin(xf*np.pi)-0.2;
ll=[pts1,pts2,pts3]
nip=[1,1]
makeMeshFromArrayList(pr1,fuselage_v,ll,nip=nip)
appendLoftedWakeRegKnownLasts(pr1,fuselage_v,1,wakeVector=np.array([30,0,0]))

fuselage_h=PanelGroup()
pts1=np.zeros((nc,3)); pts1[:,0]=xf*11-3.5;pts1[:,1]=-0.5#*np.sin(xf*np.pi)-0.2; 
pts2=np.zeros((nc,3)); pts2[:,0]=xf*11-3.5;
pts3=np.zeros((nc,3)); pts3[:,0]=xf*11-3.5;pts3[:,1]= 0.5#*np.sin(xf*np.pi)+0.2;
m1=pts1[:,0] <= 0
pts1[:,1][m1]=np.linspace(-0.1,-0.5,int(m1.sum()))
pts3[:,1][m1]=np.linspace(0.1,0.5,int(m1.sum()))
ll=[pts1,pts2,pts3]
nip=[1,1]
makeMeshFromArrayList(pr1,fuselage_h,ll,nip=nip)
appendLoftedWakeRegKnownLasts(pr1,fuselage_h,1,wakeVector=np.array([30,0,0]))



pr1.addRegion(wing_l)
pr1.addRegion(wing_r)
pr1.addRegion(htu_l)
pr1.addRegion(htu_r)
pr1.addRegion(vtu)
pr1.addRegion(fuselage_v)
pr1.addRegion(fuselage_h)
pr1.dom1.rho=1.
pr1.dom1.initSolution()
pr1.dom1.VEL=VEL


htu_l.ctrMask=np.zeros(len(htu_l.panels),dtype=bool)                                                                                                                                                    
htu_r.ctrMask=np.zeros(len(htu_r.panels),dtype=bool) 
lmHtuL=list(np.array([list(range(i*10+7,i*10+10)) for i in range(0,9)],dtype=int).flatten())
lmHtuR=list(np.array([list(range(i*10+7,i*10+10)) for i in range(1,10)],dtype=int).flatten())

htu_l.ctrMask[lmHtuL]=True  
htu_r.ctrMask[lmHtuR]=True  

htuAx1=np.array([7.35,-2.,2.5])
htuAx2=np.array([7.35,2.,2.5])

htu_l.AX_O=htuAx1
htu_l.AX=htuAx2-htuAx1
htu_r.AX_O=htuAx2
htu_r.AX=htuAx2-htuAx1


aircraft=PanelGroup();aircraft.addPanelGroup(wing_l);aircraft.addPanelGroup(wing_r);aircraft.addPanelGroup(htu_l);aircraft.addPanelGroup(htu_r);aircraft.addPanelGroup(vtu);aircraft.addPanelGroup(fuselage_h);aircraft.addPanelGroup(fuselage_v);

pr1.panelGroups={'aircraft':aircraft,'wing_l':wing_l,'wing_r':wing_r,'fuselage_h':fuselage_h,'fuselage_v':fuselage_v,'htu_l':htu_l,'htu_r':htu_r}
pr1.dom1.panelGroups=pr1.panelGroups




fcs=FlightControlSystem1()
fcs.isActive_Aileron=False                                                                                                                                                             
fcs.isActive_Rudder=False 
fcs.GRP_Rudder=[vtu]
fcs.GRP_Elevator=[htu_l,htu_r]

fd=FlightModel1()
fd.CG=cg
fd.dom1=pr1.dom1
fd.controls=fcs
pr1.dom1.COR=cg
fd.dt=1e-3

def get_FZ_MY(alpha,de):
	VEL=np.array([np.cos(alpha),0.,np.sin(alpha)])*sp
	pr1.dom1.VEL=VEL
	fcs.deflectRudders(0,de,0)
	#fcs.deflectRudders(0,0,np.deg2rad(30))
	#pr1.dom1.initSolution()
	pr1.dom1.reinitPanelsGeometry()
	pr1.dom1.updateTimeStep()
	pr1.dom1.compute()
	fz=pr1.dom1.getForce()[2]
	my=pr1.dom1.getMoment()[1]
	return np.array([fz,my])

def get_FZ_MY_trim(x,args):
	alpha=x[0]
	de=x[1]
	print(alpha,de)
	if abs(de)>np.deg2rad(25):
		de=np.sign(de)*np.deg2rad(25)
	fzReq=args[0]
	r=get_FZ_MY(alpha,de)
	return np.array([fzReq-r[0],r[1]])


print('Finding trimmed initial conditions..')
x=root(get_FZ_MY_trim,x0=np.array([0., 0.0]),args=np.array([200*9.81],),method='hybr')['x']
alpha,de=x
alphad=np.rad2deg(alpha)
ded=np.rad2deg(de)

print(f'Initial contitions are: alpha = {alphad} degrees, deltaElevator= {ded} degrees... ')
tt=                  np.array([0  ,0.1,0.4,0.6,1.2,1.5,5])
tElevator=np.deg2rad(np.array([ded,ded,-30,-30,30,ded,ded]))
tAileron=np.zeros_like(tt)
tRudder=np.zeros_like(tt)
fcs.loadTimeCmds(tt,tAileron,tElevator,tRudder)
fd.restrictedDOF.fill(False)
#fd.restrictedDOF[1]=True
#fd.restrictedDOF[3]=True
#fd.restrictedDOF[5]=True

VEL=np.array([np.cos(alpha),0.,np.sin(alpha)])*sp
pr1.dom1.VEL=VEL
	
fd.V=np.array([VEL[0],VEL[1],-VEL[2]])
fd.MASS=200.0
fd.II=np.eye(3)*1000.0

for i in tqdm(range(0,5000)):
	fcs.setDeflectionsInTime(fd.t)
	if fcs.ruddersMoved:
		#pr1.dom1.initSolution()
		pr1.dom1.reinitPanelsGeometry()
		pr1.dom1.updateTimeStep()
		
		
	fd.compute()
	fd.saveJSON('examples/ex5/JSON/time_{:d}.json'.format(i))
	pr1.plotToFile('examples/ex5/VTK/time_{:d}.vtk'.format(i),False)
	
	fd.t+=fd.dt
	
