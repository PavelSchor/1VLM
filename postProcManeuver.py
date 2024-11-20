import sys
#sys.path.append('../../')
import os.path

from collections import OrderedDict
import numpy as np
import pylab as plt
import json_tricks.np as json
import utils.csys as cs
#sys.path.append('/home/pavel/Documents/1/')

from utils.plotting import NiceGraph2D


import vtk

pathRes='analyses/mv5/man_rudder_kick_1/'
pathRes='analyses/mv5/man_elev_VA/'
pathRes='analyses/mv5/man_elev_VD/'
pathRes='analyses/mv5/man_comb_VA/'
#pathRes='analyses/mv5/man_rudder_kick_VA_beta_27deg/'

cg=np.array([2.056, 0.0, -0.102])
pHtu=np.array([6.2, 0.0 ,0.6])
pVtu=np.array([6.42, 0.0 ,1.0])

rHtu=pHtu-cg
rVtu=pVtu-cg

mHtu=10.
mVtu=10.


jds=[]
i=1
while True:
	if os.path.isfile(pathRes+'JSON/time_'+str(i)+'.json'):
		i+=1
	else:
		break

rng=range(1,i)
print('processing {:d} files'.format(i))


def makeNiceGraph(title,xlabel,ylabel):
	g1=NiceGraph2D()
	g1.xlabel=xlabel
	g1.ylabel=ylabel
	g1.title=title
	return g1


def plotCSYS2D(plt,csys,i,j):
	colors=['red','green','blue']
	for k in range(0,3):
		c=csys[k]
		plt.plot( [0,c[i]],[0,c[j]],color=colors[k])

def plotArrow2D(plt,v,i,j,scale=1.,norm=False):
	if norm:
		vv=v/np.linalg.norm(v)
	else:
		vv=v
	plt.arrow(0,0,vv[i],vv[j])
	
	
def getDictValue(jd,key,idx=None):
	if idx is None:
		return jd[key]
	else:
		return jd[key][idx]
	

def getControlInputs(jds,ix):
	x=np.zeros(len(jds))
	y=np.zeros(len(jds))
	c=0
	for d in jds:
		x[c]=getDictValue(d,'t')
		#if ix is None:
		y[c]=d['controlInputs'][ix]
		c+=1
	return x,y

def getAlphaFromDicts(jds):
	x=np.zeros(len(jds))
	y=np.zeros(len(jds))
	c=0
	for d in jds:
		x[c]=getDictValue(d,'t')
		y[c]=np.arctan2(-1*getDictValue(d,'V',2),getDictValue(d,'V',0))
		c+=1
	return x,y

def getBetaFromDicts(jds):
	x=np.zeros(len(jds))
	y=np.zeros(len(jds))
	c=0
	for d in jds:
		x[c]=getDictValue(d,'t')
		y[c]=np.arctan2(-1*getDictValue(d,'V',1),getDictValue(d,'V',0))
		c+=1
	return x,y

def getFromDicts(jds,k,ix=None):
	x=np.zeros(len(jds))
	y=np.zeros(len(jds))
	c=0
	for d in jds:
		x[c]=getDictValue(d,'t')
		if ix is None:
			y[c]=getDictValue(d,k)
		elif isinstance(ix,(int)):
			y[c]=getDictValue(d,k)[ix]
		elif isinstance(ix,tuple):
			y[c]=getDictValue(d,k)[ix[0],ix[1]]
		c+=1
	return x,y

def getData(jds,k,ix=None):
	return getFromDicts(jds,k,ix)[1]

def writeCG_POS_VTP(fname,pts):
	points=vtk.vtkPoints()
	for i in range(0,len(pts)):
		points.InsertNextPoint(pts[i])
	polydata = vtk.vtkPolyData()
	polydata.SetPoints(points)
	writer=vtk.vtkXMLPolyDataWriter()
	writer.SetFileName(fname)
	
	appendFilter=vtk.vtkAppendPolyData()
	
	appendFilter.AddInputData(polydata)
	
	#appendFilter.AddInputConnection(points.GetOutputPort())
	writer.SetInputConnection(appendFilter.GetOutputPort())
	writer.Write()

def writeCG_POS_VTP_timeSteps(fname,pts):
	for i in range(0,len(pts)):
		writeCG_POS_VTP(fname+str(i)+'.vtp',np.array([pts[i]]))
		
	
def writePVSM_array(fname,pts):
	f=open(fname,'w')
	for i in range(0,pts.shape[0]):
		for j in range(0,3):
			f.write('        <Element index="{:d}" value="{:f}"/>\n'.format(i*3+j,pts[i,j]))
	f.close()

def controlInputs2Floats(jds):
	for i in jds:
		for k in i['controlInputs']:
			i[k]=i['controlInputs'][k]
			
def dicts2Array(jds):
	nRows=len(jds)
	nCols=0
	lens=OrderedDict()
	for k in jds[0]:
		el=jds[0][k]
		if isinstance(el,float):
			lens[k]=1
			nCols+=1
		if isinstance(el,np.ndarray):
			if len(el.shape) == 1:
				lens[k]=len(el)
				nCols+=lens[k]
	res=np.zeros((nRows,nCols))
	pos=0
	sts=[]
	for k in lens:
		if lens[k]>1:
			for j in range(0,lens[k]):
				res[:,pos]=getData(jds,k,j)
				pos+=1
				sts.append(f'{k}_{j}')
		else:
			res[:,pos]=getData(jds,k)
			pos+=1
			sts.append(k)
	return res, sts

def saveAsCSV(jds,path,fname):
	res,sts=dicts2Array(jds)
	np.savetxt(f'{path}/results_{fname}.csv',res, header=str(sts))
	
	
for i in rng:
	with open(pathRes+'JSON/time_'+str(i)+'.json','r') as f:
		jds.append(json.load(f))
controlInputs2Floats(jds)
saveAsCSV(jds,pathRes,pathRes.split('/')[-2])



t,a=getAlphaFromDicts(jds)
g1=NiceGraph2D()
g1.xlabel=r'Time $  t  [s]$ '
g1.ylabel=r'Angle of attack $ \alpha  [ ^{\circ} ]$ '
g1.title='Pull up - Angle of attack '
#g1.addLine(t,np.rad2deg(a),label='Rigid')
g1.addLine(t,np.rad2deg(a),label='Computed')
g1.plot(pathRes+'t_alpha.png')



t,b=getBetaFromDicts(jds)
g1=NiceGraph2D()
g1.xlabel=r'Time $  t  [s]$ '
g1.ylabel=r'Sideslip $ beta  [ ^{\circ} ]$ '
g1.title='Sideslip '
g1.addLine(t,-np.rad2deg(b),label='Computed')
g1.plot(pathRes+'t_beta.png')




acc_s=np.zeros((len(t),3))
g_s=np.zeros((len(t),3))
t,acc_s[:,0]=getFromDicts(jds,'acceleration_s',0)
t,acc_s[:,1]=getFromDicts(jds,'acceleration_s',1)
t,acc_s[:,2]=getFromDicts(jds,'acceleration_s',2)
t,g_s[:,0]=getFromDicts(jds,'gravity_s',0)
t,g_s[:,1]=getFromDicts(jds,'gravity_s',1)
t,g_s[:,2]=getFromDicts(jds,'gravity_s',2)
g1=NiceGraph2D()
g1.xlabel=r'Time  $  t  [s]$ '
g1.ylabel=r' Accelerations_s'
g1.title=r'Accelerations_s'
g1.addLine(t,acc_s[:,0],label='$acc_s X$')
g1.addLine(t,acc_s[:,1],label='$acc_s Y$')
g1.addLine(t,acc_s[:,2],label='$acc_s Z$')

g1.addLine(t,g_s[:,0],label='$g_s X$')
g1.addLine(t,g_s[:,1],label='$g_s Y$')
g1.addLine(t,g_s[:,2],label='$g_s Z$')
g1.plot(pathRes+'t_accelerations.png')


t,omega_x=getFromDicts(jds,'angularVelocity_s',0)
t,omega_y=getFromDicts(jds,'angularVelocity_s',1)
t,omega_z=getFromDicts(jds,'angularVelocity_s',2)
g1=NiceGraph2D()
g1.xlabel=r'Time  $  t  [s]$ '
g1.ylabel=r' Angular velocity $\Omega$ [deg/s]'
g1.title=r'Angular velocity $\Omega$ [deg/s]'
g1.addLine(t,np.rad2deg(omega_x),label='$\Omega_X$')
g1.addLine(t,np.rad2deg(omega_y),label='$\Omega_Y$')
g1.addLine(t,np.rad2deg(omega_z),label='$\Omega_Z$')
g1.plot(pathRes+'t_omega.png')



t,epsilon_X=getFromDicts(jds,'angularAcceleration_s',0)
t,epsilon_Y=getFromDicts(jds,'angularAcceleration_s',1)
t,epsilon_Z=getFromDicts(jds,'angularAcceleration_s',2)
g1=NiceGraph2D()
g1.xlabel=r'Time  $  t  [s]$ '
g1.ylabel=r' Angular acceleration $\epsilon$ [deg/s^2]'
g1.title=r'Angular acceleration $\epsilon$ [deg/s^2] '
g1.addLine(t,np.rad2deg(epsilon_X),label='$\epsilon_X$')
g1.addLine(t,np.rad2deg(epsilon_Y),label='$\epsilon_Y$')
g1.addLine(t,np.rad2deg(epsilon_Z),label='$\epsilon_Z$')
g1.plot(pathRes+'t_epsilon.png')

omega_s=np.hstack([ omega_x[None].T , omega_y[None].T, omega_z[None].T])
epsilon_s=np.hstack([ epsilon_X[None].T , epsilon_Y[None].T, epsilon_Z[None].T])



inertiaForceHtu=mHtu*(g_s - (acc_s +np.cross(epsilon_s,rHtu) + +np.cross(omega_s,np.cross(omega_s,rHtu)  ))  )

inertiaForceVtu=mVtu*(g_s - (acc_s +np.cross(epsilon_s,rVtu) + +np.cross(omega_s,np.cross(omega_s,rVtu)  ))  )






t,vx=getFromDicts(jds,'V',0)
t,vz=getFromDicts(jds,'V',2)
vm=np.linalg.norm(np.vstack([vx,vz]).T,axis=1)
g1=NiceGraph2D()
g1.xlabel=r'Time  $  t  [s]$ '
g1.ylabel=r'Airspeed  $  v  [km/h]$ '
g1.title='Airspeed'
g1.addLine(t,vm*3.6,label='Computed')
g1.plot(pathRes+'t_vel.png')





t,wy=getFromDicts(jds,'OMEGA',1)
t,ey=getFromDicts(jds,'ACC_ROT',1)
#plt.plot(t,wy);plt.plot(t,ey);
#plt.show();



try:
	plt.axes(aspect='equal')
	t,x=getFromDicts(jds,'POS_EARTH',0)
	x=getFromDicts(jds,'POS_EARTH',0)[1][t>10]
	z=getFromDicts(jds,'POS_EARTH',2)[1][t>10]
	x-=x[0]
	z-=z[0]


	g1=NiceGraph2D()
	g1.xlabel=r'Position  $  x  [m]$ '
	g1.ylabel=r'Position  $  z  [m]$ '
	g1.title='Pull up - Position in Earth csys '
	g1.addLine(x,z,label='Computed')
	g1.addLine(loopsPos[:,1],loopsPos[:,2],label='Measured')
	
	plt.plot(x,z,label='Computed')
	plt.plot(loopsPos[:,1],loopsPos[:,2],label='Measured')
	plt.axes(aspect='equal')
	plt.axes(aspect='equal')

	g1.plot(pathRes+'t_pos.png')
except:
	pass

t,x=getFromDicts(jds,'POS_EARTH',0)
t,z=getFromDicts(jds,'POS_EARTH',2)
g1=NiceGraph2D()
g1.xlabel=r'Time  $  t  [s]$ '
g1.ylabel=r'Position  $  z  [m]$ '
g1.title='Pull up - Position in Earth csys '
g1.addLine(t,z,label='Computed')


g1.plot(pathRes+'t_alt.png')


#plt.axes(aspect='equal')
t,az=getFromDicts(jds,'acceleration_s',2)
t,gz=getFromDicts(jds,'gravity_s',2)
t,nz=getFromDicts(jds,'FORCE_AERO',2)

nz/=(9.8*850.)
#nz/=9.81
#nz+=1
#nz=(gz-az)/-9.81
g1=NiceGraph2D()
g1.xlabel=r'Time $  t  [s]$ '
g1.ylabel=r'Load factor $ n_z  [ - ]$ '
g1.title='Pull up - load factor $n_z$ '
g1.addLine(t,nz,label='Computed')

g1.plot(pathRes+'t_nz.png')


t,force_HTU_x=getFromDicts(jds,'force_htu',0)
t,force_HTU_y=getFromDicts(jds,'force_htu',1)
t,force_HTU_z=getFromDicts(jds,'force_htu',2)
t,force_HTU_mx=getFromDicts(jds,'force_htu',3)
t,force_HTU_my=getFromDicts(jds,'force_htu',4)
t,force_HTU_mz=getFromDicts(jds,'force_htu',5)
f_HTU=np.hstack([force_HTU_x[None].T,force_HTU_y[None].T,force_HTU_z[None].T,force_HTU_mx[None].T,force_HTU_my[None].T,force_HTU_mz[None].T,])
force_HTU=np.hstack([force_HTU_x[None].T,force_HTU_y[None].T,force_HTU_z[None].T,])
fHTU=force_HTU+inertiaForceHtu

g1=NiceGraph2D()
g1.xlabel=r'Time  $  t  [s]$ '
g1.ylabel=r' Force_HTU $F_y$'
g1.title=r'Force_HTU $F_y$ '
g1.addLine(t,force_HTU_x,label='Fx_aero')
g1.addLine(t,force_HTU_y,label='Fy_aero')
g1.addLine(t,force_HTU_z,label='Fz_aero')
g1.addLine(t,inertiaForceHtu[:,0],label='Fx_inert')
g1.addLine(t,inertiaForceHtu[:,1],label='Fy_inert')
g1.addLine(t,inertiaForceHtu[:,2],label='Fz_inert')
g1.addLine(t,np.linalg.norm(force_HTU,axis=1),label='magnitude total')
g1.plot(pathRes+'t_forces_HTU.png')


g1=NiceGraph2D()
g1.xlabel=r'Time  $  t  [s]$ '
g1.ylabel=r' Force_HTU $F_y$'
g1.title=r'Force_HTU $F_y$ '
g1.addLine(t,fHTU[:,0],label='Fx_tot')
g1.addLine(t,fHTU[:,1],label='Fy_tot')
g1.addLine(t,fHTU[:,2],label='Fz_tot')
#g1.addLine(t,np.linalg.norm(force_HTU,axis=1),label='magnitude')
g1.plot(pathRes+'t_force_total_HTU.png')



t,force_VTU_x=getFromDicts(jds,'force_vtu',0)
t,force_VTU_y=getFromDicts(jds,'force_vtu',1)
t,force_VTU_z=getFromDicts(jds,'force_vtu',2)
t,force_VTU_mx=getFromDicts(jds,'force_vtu',3)
t,force_VTU_my=getFromDicts(jds,'force_vtu',4)
t,force_VTU_mz=getFromDicts(jds,'force_vtu',5)
f_VTU=np.hstack([force_VTU_x[None].T,force_VTU_y[None].T,force_VTU_z[None].T,force_VTU_mx[None].T,force_VTU_my[None].T,force_VTU_mz[None].T,])

force_VTU=np.hstack([force_VTU_x[None].T,force_VTU_y[None].T,force_VTU_z[None].T,])
fVTU=force_VTU+inertiaForceVtu

g1=NiceGraph2D()
g1.xlabel=r'Time  $  t  [s]$ '
g1.ylabel=r' Force_VTU $F_y$'
g1.title=r'Force_VTU $F_y$ '
g1.addLine(t,force_VTU_x,label='Fx')
g1.addLine(t,force_VTU_y,label='Fy')
g1.addLine(t,force_VTU_z,label='Fz')
g1.addLine(t,np.linalg.norm(force_VTU,axis=1),label='magnitude')
g1.plot(pathRes+'t_force_VTU.png')


g1=NiceGraph2D()
g1.xlabel=r'Time  $  t  [s]$ '
g1.ylabel=r' Force_VTU $F_y$'
g1.title=r'Force_VTU $F_y$ '
g1.addLine(t,force_VTU_x,label='Fx_aero')
g1.addLine(t,force_VTU_y,label='Fy_aero')
g1.addLine(t,force_VTU_z,label='Fz_aero')
g1.addLine(t,inertiaForceVtu[:,0],label='Fx_inert')
g1.addLine(t,inertiaForceVtu[:,1],label='Fy_inert')
g1.addLine(t,inertiaForceVtu[:,2],label='Fz_inert')
g1.addLine(t,np.linalg.norm(force_VTU,axis=1),label='magnitude total')
g1.plot(pathRes+'t_forces_VTU.png')


g1=NiceGraph2D()
g1.xlabel=r'Time  $  t  [s]$ '
g1.ylabel=r' Force_VTU $F_y$'
g1.title=r'Force_VTU $F_y$ '
g1.addLine(t,fVTU[:,0],label='Fx_tot')
g1.addLine(t,fVTU[:,1],label='Fy_tot')
g1.addLine(t,fVTU[:,2],label='Fz_tot')
#g1.addLine(t,np.linalg.norm(force_VTU,axis=1),label='magnitude')
g1.plot(pathRes+'t_force_total_VTU.png')

(x,y,xlabel,ylabel,title,baseName)=(t,nz,r'Time $  t  [s]$ ',r'Load factor $ n_z  [ - ]$ ','Load factor $n_z$ ',pathRes+'data/nz/')

def arctand2(y,x):
	def arctand(y,x):
		if y>= 0 and x>=0:
			return np.rad2deg(np.arctan(y/x))
		
		if y> 0 and x<0:
			return -np.rad2deg(np.arctan(y/x))# +90
		
		if y< 0 and x<0:
			return -np.rad2deg(np.arctan(y/x))
		
		if y< 0 and x>0:
			return np.rad2deg(np.arctan(y/x))
		
	res=np.zeros_like(x)
	for i in range(0,len(x)):
		res[i]=arctand(y[i],x[i])
	return res 

t,thx=getFromDicts(jds,'CSYS_EARTH',(0,0))
thz=getFromDicts(jds,'CSYS_EARTH',(0,2))[1]
tht=(np.rad2deg(np.arctan2(thz,thx))+360)%360
tht=-arctand2(thz,thx)
#tht=np.rad2deg(np.arctan(thz/thx))


g1=NiceGraph2D()
g1.xlabel=r'Time  $  t  [s]$ '
g1.ylabel=r'Position  $  z  [m]$ '
g1.title='Pull up - Position in Earth csys '
g1.addLine(t,tht,label='Computed')

plt.axes(aspect='equal')

g1.plot(pathRes+'t_pitch.png')



t,thy=getFromDicts(jds,'CSYS_EARTH',(2,1))
thz=getFromDicts(jds,'CSYS_EARTH',(2,2))[1]
phi=-np.rad2deg(np.arctan2(thy,thz))
#tht=np.rad2deg(np.arctan(thz/thx))
#tht[np.where((tht > 90) & (tht < 180), True, False)]= 180 -tht[np.where((tht > 90) & (tht < 180), True, False)]
#tht[np.where((tht >= 180) & (tht < 270), True, False)]= 180+tht[np.where((tht >= 180) & (tht < 270), True, False)]
t,phi=getFromDicts(jds,'RPOS_EARTH',0)
phi=np.rad2deg(phi)
t,y=getFromDicts(jds,'POS_EARTH',1)
t,z=getFromDicts(jds,'POS_EARTH',2)
g1=NiceGraph2D()
g1.xlabel=r'Time  $  t  [s]$ '
g1.ylabel=r'P  $  \Phi  [deg]$ '
g1.title='Bank angle '
g1.addLine(t,-phi,label='Computed')
plt.axes(aspect='equal')

g1.plot(pathRes+'t_roll.png')







t,dAileronL=getControlInputs(jds,'dAileron')
t,dElevator=getControlInputs(jds,'dElevator')
t,dRudder=getControlInputs(jds,'dRudder')
g1=NiceGraph2D()
g1.xlabel=r'time  $  t  [s]$ '
g1.ylabel=r'deflections  []$ '
g1.title='deflections'
g1.addLine(t,dAileronL,label='dAileronL')
g1.addLine(t,dElevator,label='dElevator')
g1.addLine(t,dRudder,label='dRudder')
g1.plot(pathRes+'t_controls.png')





#for i in range(0,x.shape[0]):
	#plt.clf()
	
	#plt.xlabel=xlabel#r'Time $  t  [s]$ '
	#plt.ylabel=ylabel#r'Load factor $ n_z  [ - ]$ '
	#plt.title=title#'Pull up - load factor $n_z$ '
	#plt.xlim(0,16)
	#plt.ylim(-1,9)
	#plt.plot(x[0:i],y[0:i])
	#plt.savefig(baseName+'t_{:d}.png'.format(i))


#(x,y,xlabel,ylabel,title,baseName)=(t,np.rad2deg(a),r'Time $  t  [s]$ ',r'Load factor $ n_z  [ - ]$ ','Angle of attack $n_z$ ',pathRes+'data/alpha/')

#for i in range(0,x.shape[0]):
	#plt.clf()
	
	#plt.xlabel=xlabel#r'Time $  t  [s]$ '
	#plt.ylabel=ylabel#r'Load factor $ n_z  [ - ]$ '
	#plt.title=title#'Pull up - load factor $n_z$ '
	#plt.xlim(0,16)
	#plt.ylim(-15,20)
	#plt.plot(x[0:i],y[0:i])
	#plt.savefig(baseName+'t_{:d}.png'.format(i))




t,my=getFromDicts(jds,'FORCE_AERO',4)
my[abs(my)>1000]=0.
g1=NiceGraph2D()
g1.xlabel=r'Time $  t  [s]$ '
g1.ylabel=r'Pitchin moment $ m_y  [Nm ]$ '
g1.title='Pull up - momeny $m_y$ '
g1.addLine(t,my,label='Rigid')
g1.plot(pathRes+'t_my.png')

#t,dEl=getFromDicts(jds,'elevator')
#g1=NiceGraph2D()
#g1.xlabel=r'Time $  t  [s]$ '
#g1.ylabel=r'Elevator deflection $ \delta_e  [ ^{\circ} ]$ '
#g1.title='Pull up - Elevator deflection $ \delta_e $ '
#g1.addLine(t,dEl,label='Rigid')
#g1.plot(pathRes+'t_elev.png')

######################## BACK
pos1=np.array([-20,0,7])
focus1=np.array([15,0,-2])
upview1=np.array([0,0,1])
CG=np.zeros((len(jds),3))
fcp=np.zeros((len(jds),3))
pos=np.zeros((len(jds),3))
vup=np.zeros((len(jds),3))

for i in range(0,len(jds)):
	CG_POS_G=jds[i]['POS_EARTH']
	cgi=jds[i]['CG']
	csE=jds[i]['CSYS_EARTH']
	csB=jds[i]['CSYS_BODY']
	pos[i]=cs.transformVector(pos1,csE,csB)+CG_POS_G-cgi
	fcp[i]=cs.transformVector(focus1,csE,csB)+CG_POS_G-cgi
	vup[i]=cs.transformVector(upview1,csE,csB)
	

np.savetxt(pathRes+'camFocalPoint.txt',fcp)
np.savetxt(pathRes+'camPosition.txt',pos)
np.savetxt(pathRes+'camViewUp.txt',vup)


######################## BACK
pos1=np.array([-20,0,7])
focus1=np.array([15,0,-2])
upview1=np.array([0,0,1])
CG=np.zeros((len(jds),3))
fcp=np.zeros((len(jds),3))
pos=np.zeros((len(jds),3))
vup=np.zeros((len(jds),3))

for i in range(0,len(jds)):
	CG_POS_G=jds[i]['POS_EARTH']
	cgi=jds[i]['CG']
	csE=jds[i]['CSYS_EARTH']
	csB=jds[i]['CSYS_BODY']
	pos[i]=cs.transformVector(pos1,csE,csB)+CG_POS_G-cgi
	fcp[i]=cs.transformVector(focus1,csE,csB)+CG_POS_G-cgi
	vup[i]=cs.transformVector(upview1,csE,csB)
	

np.savetxt(pathRes+'camFocalPoint.txt',fcp)
np.savetxt(pathRes+'camPosition.txt',pos)
np.savetxt(pathRes+'camViewUp.txt',vup)


######################## LEFT TOP
pos1=np.array([6,-14,6])
focus1=np.array([-17,13,-4])
upview1=np.array([0,0,1])
CG=np.zeros((len(jds),3))
fcp=np.zeros((len(jds),3))
pos=np.zeros((len(jds),3))
vup=np.zeros((len(jds),3))

for i in range(0,len(jds)):
	CG_POS_G=jds[i]['POS_EARTH']
	cgi=jds[i]['CG']
	csE=jds[i]['CSYS_EARTH']
	csB=jds[i]['CSYS_BODY']
	pos[i]=cs.transformVector(pos1,csE,csB)+CG_POS_G-cgi
	fcp[i]=cs.transformVector(focus1,csE,csB)+CG_POS_G-cgi
	vup[i]=cs.transformVector(upview1,csE,csB)
	

np.savetxt(pathRes+'camFocalPoint.txt',fcp)
np.savetxt(pathRes+'camPosition.txt',pos)
np.savetxt(pathRes+'camViewUp.txt',vup)



######################## LEFT TOP FIXED Z
pos1=np.array([6,-20,6])
focus1=np.array([-17,13,-4])
upview1=np.array([0,0,1])
CG=np.zeros((len(jds),3))
fcp=np.zeros((len(jds),3))
pos=np.zeros((len(jds),3))
vup=np.zeros((len(jds),3))

for i in range(0,len(jds)):
	CG_POS_G=jds[i]['POS_EARTH']
	cgi=jds[i]['CG']
	csE=jds[i]['CSYS_EARTH']
	csB=jds[i]['CSYS_BODY']
	pos[i]=cs.transformVector(pos1,csE,csB)+CG_POS_G-cgi
	fcp[i]=cs.transformVector(focus1,csE,csB)+CG_POS_G-cgi
	vup[i]=cs.transformVector(upview1,csE,csB)
	
vup[:,0].fill(0)
vup[:,1].fill(0)
vup[:,2].fill(1)

np.savetxt(pathRes+'camFocalPoint.txt',fcp)
np.savetxt(pathRes+'camPosition.txt',pos)
np.savetxt(pathRes+'camViewUp.txt',vup)
