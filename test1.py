from utils.mesh import *
from vlm.vlmMethod1 import *
import numpy as np
np.setbufsize(8192*2)
import pylab as plt
import time


chord=1.
span=6.
alpha=np.deg2rad(10)
sp=35.
q=0.5*1*sp**2
VEL=np.array([np.cos(alpha),0,np.sin(alpha)])*sp
pr1=VlmProblem()

############wing=PanelGroup()
############makeSpanwiseLoftedPlate1Reg(pr1,wing,50,10,np.array([0,-span/2,0]),np.array([0,span/2,0]),np.array([chord,span/2,0]),np.array([chord,-span/2,0]))
############appendLoftedWakeReg(pr1,wing,10,wakeVector=np.array([100,0,0]))

#############zw=10
#############wing2=PanelGroup()
#############makeSpanwiseLoftedPlate1Reg(pr1,wing2,100,10,np.array([0,-span/2,zw]),np.array([0,span/2,zw]),np.array([chord,span/2,zw]),np.array([chord,-span/2,zw]))
#############appendLoftedWakeReg(pr1,wing2,100,wakeVector=np.array([100,0,0]))

wing=PanelGroup()
makeSpanwiseLoftedPlate1Reg(pr1,wing,40,10,np.array([0,-span/2,0]),np.array([0,0,0]),np.array([chord,0,0]),np.array([chord,-span/2,0]))
appendLoftedWakeReg(pr1,wing,10,wakeVector=np.array([100,0,0]))

zw=0.0
wing2=PanelGroup()
makeSpanwiseLoftedPlate1Reg(pr1,wing2,40,10,np.array([0,0,zw]),np.array([0,span/2,zw]),np.array([chord,span/2,zw]),np.array([chord,0,zw]))
appendLoftedWakeReg(pr1,wing2,10,wakeVector=np.array([100,0,0]))


pr1.addRegion(wing)
pr1.addRegion(wing2)

#self.computeWakeContribution()
chords=np.ones(pr1.dom1.nPan)*chord
VTOT=np.tile(VEL,(pr1.dom1.nPan,1))
pr1.dom1.rho=1.
pr1.dom1.initSolution()
pr1.dom1.VEL=VEL

#%time 
pr1.dom1.compute()

#pr1.addRegion(wake)



self=pr1.dom1


#for i in self.panels:
	#self.panels[i].G=self.GAMA[i]
	#self.panels[i].pressureCoef=self.GAMA[i]#*1.225*sp
#region=self.regions[0]
#m=region.nSpanwise 
#n=region.nChordwise
#for i in range(0,m):
	#for j in range(1,n):
		#p=region.panels[i*n+j]
		#pp=region.panels[i*n+j-1]
		#p.pressureCoef=p.G-pp.G
xyz=np.zeros((len(self.panels),3))

cf=np.zeros(len(self.panels))
for i in range(0, len(self.panels)):
	xyz[i]=self.panels[i].getCpoint()
	cf[i]=self.panels[i].pressureCoef
plt.plot(xyz[:,0],cf);plt.show();
plt.plot(xyz[:,1],cf);plt.show();

#self.FORCE.sum()/q/chord/span
CF=self.FORCE.sum(axis=0)/q/chord/span
CL=-CF[0]*np.sin(alpha)+CF[2]*np.cos(alpha)



alphaRng=np.deg2rad(np.linspace(-0,20,21))


CL=np.zeros_like(alphaRng)
CD=np.zeros_like(alphaRng)
CM=np.zeros_like(alphaRng)

start = time.time()
for i in range(0,len(alphaRng)):
	alpha=alphaRng[i]
	VEL=np.array([np.cos(alpha),0.0,np.sin(alpha)])*sp
	pr1.dom1.VEL=VEL


	#%time pr1.dom1.compute()
	#%time pr1.dom1.computeViscous()
	#%timeit pr1.dom1.compute()


	##%timeit 
	#pr1.dom1.computeViscous(10)
	
	#pr1.dom1.compute()
	#%time 
	pr1.dom1.updateTimeStep()
	pr1.dom1.compute()
	#print abs(pr1.dom1.RESIDUAL).max()
	CF=self.FORCE.sum(axis=0)/q/chord/span
	CL[i]=-CF[0]*np.sin(alpha)+CF[2]*np.cos(alpha)
	CD[i]= CF[0]*np.cos(alpha)+CF[2]*np.sin(alpha)
	CM[i]=self.MOMENT[1]/q/chord/span/chord
	pr1.plotToFile('figs/wing_plain{:d}.vtp'.format(i),True)
	fw=wing.getForceAtRefPT()
	
	print(np.rad2deg(alpha), ': ', fw[3]/fw[2])
#end = time.time()
#print("Elapsed (after compilation) = %s" % (end - start))


#region=self.regions[0]
#ctrMask=np.zeros(region.nPanels,dtype=np.bool)
#for i in region.panels:
	#p=region.panels[i]
	#if np.all(p.pts[:,0]>=0.7):
		#ctrMask[i]=True
##	def deflectControlSurface(self,ctrMask,ax,orig,d):

#region.deflectControlSurface(ctrMask=ctrMask,ax=np.array([0,1,0]),orig=np.array([0.7,0,0]),d=np.deg2rad(30))

#self.initSolution()
#self.reinitPanelsGeometry()
#self.compute()
#pr1.plotToFile('a.vtp',True)



#CL_F=np.zeros_like(alphaRng)
#CD_F=np.zeros_like(alphaRng)
#CM_F=np.zeros_like(alphaRng)

#start = time.time()
#for i in range(0,len(alphaRng)):
	#alpha=alphaRng[i]
	#VEL=np.array([np.cos(alpha),0.0,np.sin(alpha)])*sp
	#pr1.dom1.VEL=VEL


	##%time pr1.dom1.compute()
	##%time pr1.dom1.computeViscous()
	##%timeit pr1.dom1.compute()


	###%timeit 
	##pr1.dom1.computeViscous(10)
	
	##pr1.dom1.compute()
	##%time 
	#pr1.dom1.updateTimeStep()
	#pr1.dom1.compute()
	##print abs(pr1.dom1.RESIDUAL).max()
	#CF=self.FORCE.sum(axis=0)/q/chord/span
	#CL_F[i]=-CF[0]*np.sin(alpha)+CF[2]*np.cos(alpha)
	#CD_F[i]= CF[0]*np.cos(alpha)+CF[2]*np.sin(alpha)
	#CM_F[i]=self.MOMENT[1]/q/chord/span/chord
	#pr1.plotToFile('figs/wing_flap{:d}.vtp'.format(i),True)
#end = time.time()
#print("Elapsed (after compilation) = %s" % (end - start))
