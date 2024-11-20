import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
sns.set_style("whitegrid")
import json_tricks as json
import matplotlib.lines as lines
from utils.beamCalc import BeamCrossSection
import time
from am_analysis.amFunctions import *
from flightLoads.pointLoads import *
from inertia.inertiaModel1 import InertiaModel
from utils.b3Vect import B3vect 
from utils.bCad import BCad
from resampleFoil import resample
from utils.bCad import *

bc=BCad()


class GraphLine2D(object):
	pass

class NiceGraph2D(object):
	def __init__(self):
		self.lines=[]
		self.xlabel=''
		self.ylabel=''
		self.title=''
		self.legendLoc=2
		self.markers = list(lines.Line2D.filled_markers)
		self.markers.sort()
		self.n=0
		pass
	
	def addLine(self,x,y=None,label='',linewidth=1):
		l=GraphLine2D()
		l.ID=self.n
		self.n+=1
		l.x=x
		if not y is None:
			l.y=y
		else:
			l.y=x[:,1]
		l.label=label
		l.linewidth=linewidth
		self.lines.append(l)
	
	def plot(self,fname=None):
		plt.clf()
		
	
		fig = plt.figure()
		ax = plt.subplot(111)


		for l in self.lines:
			
			ax.plot(l.x,l.y,label=l.label,linewidth=l.linewidth,marker = self.markers[l.ID%len(self.markers)],markersize=5)


		plt.grid(True)
		plt.xlabel(self.xlabel)
		plt.ylabel(self.ylabel)
		plt.title(self.title)
		# Shrink current axis by 20%
		box = ax.get_position()
		ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

		# Put a legend to the right of the current axis
		ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),ncol=2)#,mode='expand')
		plt.setp(plt.gca().get_legend().get_texts(), fontsize='6')
		if not fname is None:
			plt.savefig(fname, dpi=400)
		else:
			plt.show()
		plt.close()	


def writeTableToHTMLStream(f,table,head=None):
	if head is None:
		head = ['' for i in range(0,table.shape[1])]
	
	f.write('<table style="width:80%">')
	f.write('<tr>\n')
	for i in head:
		f.write('<th>{:s}</th>'.format(str(i)))
	f.write('</tr>\n')
	for i in range(0,table.shape[0]):
		f.write('<tr>\n')
		for j in range(0,table.shape[1]):
			f.write('<th>{:3.1f}</th>'.format(table[i][j]))
		f.write('</tr>\n')
	f.write('</table>\n')


def postProcWingLoad(path,cases,folderName,sectionName):
	ic=1
	ta=NiceGraph2D()
	ta.xlabel='Span y [m]'
	ta.ylabel='Force [N]'
	ta.title='Fy axial forces'

	tn=NiceGraph2D()
	tn.xlabel='Span y [m]'
	tn.ylabel='Force [N]'
	tn.title='Fz normal forces'
	
	tt=NiceGraph2D()
	tt.xlabel='Span y [m]'
	tt.ylabel='Force [N]'
	tt.title='Fx tangential forces'
	
	mn=NiceGraph2D()
	mn.xlabel='Span y [m]'
	mn.ylabel='Moment [Nm]'
	mn.title='Mx normal bending moments'

	mt=NiceGraph2D()
	mt.xlabel='Span y [m]'
	mt.ylabel='Moment [Nm]'
	mt.title='Mz tangential bending moments'
	
	mk=NiceGraph2D()
	mk.xlabel='Span y [m]'
	mk.ylabel='Moment [Nm]'
	mk.title='My torsion moments'
	
	
	for cid in cases:
		fname=path+str(cid)+'/definition.json'
		with open(fname,'r') as data_file:
			d = json.load(data_file)
		pathRes=d['path']
		wl=np.loadtxt(pathRes+folderName+sectionName+'_LoadTotal.txt')
		wla=np.loadtxt(pathRes+folderName+sectionName+'_LoadAero.txt')
		wlm=np.loadtxt(pathRes+folderName+sectionName+'_LoadInertia.txt')
		g1=NiceGraph2D()
		g1.xlabel='Span y [m]'
		g1.ylabel='Force [N]'
		g1.title=d['description']+' forces'
		
		g2=NiceGraph2D()
		g2.xlabel='Span y [m]'
		g2.ylabel='Moment [Nm]'
		g2.title=d['description']+' moments'
		
		g1.addLine(wl[:,ic],wl[:,3],label='Fx')
		g1.addLine(wl[:,ic],wl[:,4],label='Fy')
		g1.addLine(wl[:,ic],wl[:,5],label='Fz')
		
		g2.addLine(wl[:,ic],wl[:,6],label='Mx')
		g2.addLine(wl[:,ic],wl[:,7],label='My')
		g2.addLine(wl[:,ic],wl[:,8],label='Mz')
		#g1.plot(pathRes+'wing/wingLoad_forces.png')
		#g2.plot(pathRes+'wing/wingLoad_moments.png')

		tt.addLine(wl[:,ic],wl[:,3],label=''+str(d['ID']))
		ta.addLine(wl[:,ic],wl[:,4],label=''+str(d['ID']))
		tn.addLine(wl[:,ic],wl[:,5],label=''+str(d['ID']))
		
		mn.addLine(wl[:,ic],wl[:,6],label=''+str(d['ID']))
		mk.addLine(wl[:,ic],wl[:,7],label=''+str(d['ID']))
		mt.addLine(wl[:,ic],wl[:,8],label=''+str(d['ID']))

	tt.plot(path+'loadEnvelope/'+folderName+sectionName+'_Load_fx.png')
	ta.plot(path+'loadEnvelope/'+folderName+sectionName+'_Load_fy.png')
	tn.plot(path+'loadEnvelope/'+folderName+sectionName+'_Load_fz.png')
	mn.plot(path+'loadEnvelope/'+folderName+sectionName+'_Load_mx.png')
	mk.plot(path+'loadEnvelope/'+folderName+sectionName+'_Load_my.png')
	mt.plot(path+'loadEnvelope/'+folderName+sectionName+'_Load_mz.png')
	#postProcWingLoad('am_analysis/loads/cases/',[i for i in range(0,8)])



def postProcWingLoad2(path,cases,folderName,sectionName):
	ic=1
	for cid in cases:
		ta=NiceGraph2D()
		ta.xlabel='Span y [m]'
		ta.ylabel='Force [N]'
		ta.title='Fy axial forces'

		tn=NiceGraph2D()
		tn.xlabel='Span y [m]'
		tn.ylabel='Force [N]'
		tn.title='Fz normal forces'
		
		tt=NiceGraph2D()
		tt.xlabel='Span y [m]'
		tt.ylabel='Force [N]'
		tt.title='Fx tangential forces'
		
		mn=NiceGraph2D()
		mn.xlabel='Span y [m]'
		mn.ylabel='Moment [Nm]'
		mn.title='Mx normal bending moments'

		mt=NiceGraph2D()
		mt.xlabel='Span y [m]'
		mt.ylabel='Moment [Nm]'
		mt.title='Mz tangential bending moments'
		
		mk=NiceGraph2D()
		mk.xlabel='Span y [m]'
		mk.ylabel='Moment [Nm]'
		mk.title='My torsion moments'
	
		fname=path+str(cid)+'/definition.json'
		with open(fname,'r') as data_file:
			d = json.load(data_file)
		pathRes=d['path']

		writeCaseReport(d)

		wl=np.loadtxt(pathRes+folderName+sectionName+'_LoadTotal.txt')
		wla=np.loadtxt(pathRes+folderName+sectionName+'_LoadAero.txt')
		wlm=np.loadtxt(pathRes+folderName+sectionName+'_LoadInertia.txt')

		tt.addLine(wl[:,ic],wl[:,3],label='Fx tangential')
		tt.addLine(wla[:,ic],wla[:,3],label='Fx aero')
		tt.addLine(wla[:,ic],wlm[:,3],label='Fx inert.')
		
		ta.addLine(wl[:,ic],wl[:,4],label='Fy axial')
		ta.addLine(wla[:,ic],wla[:,4],label='Fy aero')
		ta.addLine(wlm[:,ic],wlm[:,4],label='Fy inert.')
		
		tn.addLine(wl[:,ic],wl[:,5],label='Fz Nornal')
		tn.addLine(wla[:,ic],wla[:,5],label='Fz aero')
		tn.addLine(wlm[:,ic],wlm[:,5],label='Fz inert.')
		
		mn.addLine(wl[:,ic],wl[:,6],label='Mx')
		mn.addLine(wla[:,ic],wla[:,6],label='Mx aero')
		mn.addLine(wlm[:,ic],wlm[:,6],label='Mx inert.')
		
		mk.addLine(wl[:,ic],wl[:,7],label='My')
		mk.addLine(wla[:,ic],wla[:,7],label='My aero')
		mk.addLine(wlm[:,ic],wlm[:,7],label='My inert.')
		
		mt.addLine(wl[:,ic],wl[:,8],label='Mz')
		mt.addLine(wla[:,ic],wla[:,8],label='Mz aero')
		mt.addLine(wlm[:,ic],wlm[:,8],label='Mz inert.')

		tt.plot(pathRes+folderName+sectionName+'_Load_fx.png')
		ta.plot(pathRes+folderName+sectionName+'_Load_fy.png')
		tn.plot(pathRes+folderName+sectionName+'_Load_fz.png')
		mn.plot(pathRes+folderName+sectionName+'_Load_mx.png')
		mk.plot(pathRes+folderName+sectionName+'_Load_my.png')
		mt.plot(pathRes+folderName+sectionName+'_Load_mz.png')
	#postProcWingLoad('am_analysis/loads/cases/',[i for i in range(0,8)])



def writeCaseReport(d,prefix='wing_left_'):
	with open(d['path']+'results/integralForces.json') as data_file:
		intForces = json.load(data_file)
		
	f=open(d['path']+'results/report.html','w')
	f.writelines('''<!DOCTYPE html>
<html>
<head>''')

	f.write('<title>Load Case'+str(d['ID'])+'</title>')

	f.writelines('''<style>
table, th, td {
    border: 1px solid black;
    border-collapse: collapse;
}
th, td {x
    padding: 5px;
    text-align: left;
}
</style>
</head>
<body>

\n''')

	f.write('<h1>Load Case '+str(d['ID'])+' report</h1>\n')
	f.write('<br> (c) 2017 Pavel Schor, pavel.schor@email.cz <br> <br><br>\n')
	f.write('THIS DOCUMENT WAS GENERATED ON : '+time.strftime('%A %d-%m-%Y - %R '))
	f.write('<h2> Case summary table </h2>')
	f.write('<table style="width:80%">')
	f.write('<tr><th>{:s}</th><th>{:d}</th></tr>\n'''.format('ID:  ',d['ID']))
	complianceCodes=str([ ss.encode('ascii','ignore') for ss in d['complianceCodes']]).replace('[','').replace(']','').replace(',','\n').replace("'","")

	f.write('<tr><th>{:s}</th><th>{:s}</th></tr>\n'''.format('Airworhiness compliace:',complianceCodes))
	f.write('<tr><th>{:s}</th><th>{:s}</th></tr>\n'''.format('Description:  ',d['description']))
	f.write('<tr><th>{:s}</th><th>{:3.1f}</th></tr>\n'.format('EAS [m/s]:  ',d['EAS']))
	f.write('<tr><th>{:s}</th><th>{:3.1f}</th></tr>\n'.format('TAS [m/s]:  ',d['TAS']))
	f.write('<tr><th>{:s}</th><th>{:3.1f}</th></tr>\n'.format('H [m]:  ',d['H']))
	f.write('<tr><th>{:s}</th><th>{:3.2f}</th></tr>\n'.format('n [-]: ',d['n']))	
	f.write('<tr><th>{:s}</th><th>{:4.2f}</th></tr>\n'.format('m [kg]: ',d['m']))
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('rho [kg*m^-3]:  ',d['rho']))
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('alpha [deg]:  ',np.rad2deg(d['alpha'])))
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('beta [deg]:  ',np.rad2deg(d['beta'])))
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('Xcg [m]:  ',d['CG'][0]))
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('Ycg [m]:  ',d['CG'][1]))
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('Zcg [m]:  ',d['CG'][2]))
	f.write('<tr><th colspan="2">   </th></tr>\n')


	f.write('<tr><th colspan="2">Linear acceleration at CG:</th></tr>\n')
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('aX [m/s^2]:  ',d['linearAcceleration'][0]))
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('aY [m/s^2]:  ',d['linearAcceleration'][1]))
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('aZ [m/s^2]:  ',d['linearAcceleration'][2]))
	f.write('<tr><th colspan="2">   </th></tr>\n')

	f.write('<tr><th colspan="2">Angular acceleration at CG:</th></tr>\n')
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('EpX [rad/s^2]:  ',d['angularAcceleration'][0]))
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('EpY [rad/s^2]:  ',d['angularAcceleration'][1]))
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('EpZ [rad/s^2]:  ',d['angularAcceleration'][2]))
	f.write('<tr><th colspan="2">   </th></tr>\n')

	f.write('<tr><th colspan="2">Angular velocity at CG:</th></tr>\n')
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('wX [rad/s]:  ',d['angularVelocity'][0]))
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('wY [rad/s]:  ',d['angularVelocity'][1]))
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('wZ [rad/s]:  ',d['angularVelocity'][2]))
	f.write('<tr><th colspan="2">   </th></tr>\n')
	
	
	f.write('<tr><th colspan="2">Aerodynamic Force FA:</th></tr>\n')
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('FXaero [N]:  ',d['aerodynamicForce'][0]))
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('FYaero [N]:  ',d['aerodynamicForce'][1]))
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('FZaero [N]:  ',d['aerodynamicForce'][2]))
	f.write('<tr><th colspan="2">   </th></tr>\n')
	
	f.write('<tr><th colspan="2">Aerodynamic Force WING FW:</th></tr>\n')
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('FWXaero [N]:  ',intForces['wingForce'][0]))
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('FWYaero [N]:  ',intForces['wingForce'][1]))
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('FWZaero [N]:  ',intForces['wingForce'][2]))

	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('MWXaero [Nm]:  ',intForces['wingForce'][3]))
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('MWYaero [Nm]:  ',intForces['wingForce'][4]))
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('MWZaero [Nm]:  ',intForces['wingForce'][5]))
	f.write('<tr><th colspan="2">   </th></tr>\n')
	
	f.write('<tr><th colspan="2">Aerodynamic Force FUSELAGE FF:</th></tr>\n')
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('FFXaero [N]:  ',intForces['fuselageForce'][0]))
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('FFYaero [N]:  ',intForces['fuselageForce'][1]))
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('FFZaero [N]:  ',intForces['fuselageForce'][2]))

	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('MFXaero [Nm]:  ',intForces['fuselageForce'][3]))
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('MFYaero [Nm]:  ',intForces['fuselageForce'][4]))
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('MFZaero [Nm]:  ',intForces['fuselageForce'][5]))
	f.write('<tr><th colspan="2">   </th></tr>\n')
	
	f.write('<tr><th colspan="2">Tail balancing Force FT:</th></tr>\n')
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('FXt [N]:  ',d['tailForce'][0]))
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('FYt [N]:  ',d['tailForce'][1]))
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('FZt [N]:  ',d['tailForce'][2]))
	f.write('<tr><th colspan="2">   </th></tr>\n')
	
	
	f.write('<tr><th colspan="2">Thrust force FTH:</th></tr>\n')
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('FXth [N]:  ',d['thrustForce'][0]))
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('FYth [N]:  ',d['thrustForce'][1]))
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('FZth [N]:  ',d['thrustForce'][2]))
	f.write('<tr><th colspan="2">   </th></tr>\n')
	
	f.write('<tr><th colspan="2">Inertia Force FI:</th></tr>\n')
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('FXin [N]:  ',d['inertiaForce'][0]))
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('FYin [N]:  ',d['inertiaForce'][1]))
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('FZin [N]:  ',d['inertiaForce'][2]))
	f.write('<tr><th colspan="2">   </th></tr>\n')
	
	f.write('<tr><th colspan="2">Viscous drag Force Fvd:</th></tr>\n')
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('FXvd [N]:  ',d['viscousDragForce'][0]))
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('FYvd [N]:  ',d['viscousDragForce'][1]))
	f.write('<tr><th>{:s}</th><th>{:4.3f}</th></tr>\n'.format('FZvd [N]:  ',d['viscousDragForce'][2]))
	f.write('</table>\n')

	
	f.write('<h2>Forces overview</h2>\n')
	f.write('<br>')
	f.write('<img src="../definition.png "alt="LOAD" width="1024px" height="auto">\n')# style="width:254px;height:86px;">'''
	f.write('<br>')
	f.write('<h2>Wing load</h2>')
	f.write('<br>')
	f.write('<h3>Wing local forces TT,TA,TN</h3>\n')
	f.write('<p>Forces and moments acting at main spar of the wing transformed to local cocrdinate system of the wing</p>')
	f.write('<p>Fx - tangential force; Fy - axial force; Fz - normal force;</p>')
	f.write('<p>Coordinates of main spar [x,y,z] are given in tables bellow</p>')
	
	f.write('<br>')
	f.write('<img src="../wing/'+prefix+'Load_fx.png" alt="TT"  width="1024px" height="auto">\n')# style="width:254px;height:86px;">'''
	f.write('<img src="../wing/'+prefix+'Load_fy.png" alt="TA"  width="1024px" height="auto">\n')# style="width:254px;height:86px;">'''
	f.write('<img src="../wing/'+prefix+'Load_fz.png" alt="TN"  width="1024px" height="auto">\n')# style="width:254px;height:86px;">'''
	f.write('<br>')
	f.write('<h3>Wing local moments MN,MK,MT\n</h3>')
	f.write('<p>Moments acting at main spar of the wing transformed to local cocrdinate system of the wing</p>')
	f.write('<p>Mx -normal bending moment; My - torsion moment; Mz - tangential bending moment;</p>')
	f.write('<p>Coordinates of main spar [x,y,z] are given in tables bellow</p>')
	
	f.write('<br>')
	f.write('<img src="../wing/'+prefix+'Load_mx.png" alt="MN"  width="1024px" height="auto">\n')# style="width:254px;height:86px;">'''
	f.write('<img src="../wing/'+prefix+'Load_my.png" alt="MK"  width="1024px" height="auto">\n')# style="width:254px;height:86px;">'''
	f.write('<img src="../wing/'+prefix+'Load_mz.png" alt="MT"  width="1024px" height="auto">\n')# style="width:254px;height:86px;">'''
	f.write('<br>')
	
	pathRes=d['path']
	f.write('<h3>Wing total forces</h3>\n')
	f.write('<p>Total forces and moments acting at main spar of the wing transformed to local cocrdinate system of the wing</p>')
	writeTableToHTMLStream(f,np.loadtxt(pathRes+'wing/'+prefix+'LoadTotal.txt'),['x','y','z','Fx','Fy','Fz','Mx','My','Mz'])
	f.write('<h3>Aerodynamic forces</h3>\n')
	f.write('<p>Aerodynamic forces and moments acting at main spar of the wing transformed to local cocrdinate system of the wing</p>')
	writeTableToHTMLStream(f,np.loadtxt(pathRes+'wing/'+prefix+'LoadAero.txt'),['x','y','z','Fx','Fy','Fz','Mx','My','Mz'])
	f.write('<h3>Inertia forces</h3>\n')
	f.write('<p>Inertia forces and moments acting at main spar of the wing transformed to local cocrdinate system of the wing</p>')
	writeTableToHTMLStream(f,np.loadtxt(pathRes+'wing/'+prefix+'LoadInertia.txt'),['x','y','z','Fx','Fy','Fz','Mx','My','Mz'])
	
	f.write('</body>\n')
	f.write('</html>\n')
	f.close()

def makeWingLoadEnvelope(path,cases,folderName,sectionName):
	cid=cases[0]
	fname=path+str(cid)+'/definition.json'
	with open(fname,'r') as data_file:
		d = json.load(data_file)
	pathRes=d['path']
	wlEnvPlus=np.loadtxt(pathRes+folderName+sectionName+'_LoadTotal.txt')
	wlEnvMinus=np.loadtxt(pathRes+folderName+sectionName+'_LoadTotal.txt')
	wlEnvAbs=np.loadtxt(pathRes+folderName+sectionName+'_LoadTotal.txt')
	wlEnvID=np.loadtxt(pathRes+folderName+sectionName+'_LoadTotal.txt')
	wlEnvID[:,3:10]=0
	wlEnvAbsXYZ=np.loadtxt(pathRes+folderName+sectionName+'_LoadTotalXYZ.txt')
	wlEnvIDXYZ=np.loadtxt(pathRes+folderName+sectionName+'_LoadTotalXYZ.txt')
	wlEnvIDXYZ[:,3:10]=0
	for cid in cases:
		fname=path+str(cid)+'/definition.json'
		with open(fname,'r') as data_file:
			d = json.load(data_file)
		pathRes=d['path']
		wl=np.loadtxt(pathRes+folderName+sectionName+'_LoadTotal.txt')
		wlXYZ=np.loadtxt(pathRes+folderName+sectionName+'_LoadTotalXYZ.txt')
		for i in range(0,wl.shape[0]):
			for j in range(3,9):
				if abs(wl[i,j])>abs(wlEnvAbs[i,j]):
					wlEnvAbs[i,j]=wl[i,j]
					wlEnvID[i,j]=cid
				if abs(wlXYZ[i,j])>abs(wlEnvAbsXYZ[i,j]):
					wlEnvAbsXYZ[i,j]=wlXYZ[i,j]
					wlEnvIDXYZ[i,j]=cid
				if abs(wl[i,j])>abs(wlEnvPlus[i,j]):
					wlEnvPlus[i,j]=wl[i,j]
				if abs(wl[i,j])>abs(wlEnvMinus[i,j]):
					wlEnvMinus[i,j]=wl[i,j]
	np.savetxt(path+'loadEnvelope/'+folderName+sectionName+'LoadEnvelopeAbs.txt',wlEnvAbs)
	np.savetxt(path+'loadEnvelope/'+folderName+sectionName+'LoadEnvelopePlus.txt',wlEnvPlus)
	np.savetxt(path+'loadEnvelope/'+folderName+sectionName+'LoadEnvelopeMinus.txt',wlEnvMinus)
	np.savetxt(path+'loadEnvelope/'+folderName+sectionName+'LoadEnvelopeAbsXYZ.txt',wlEnvAbsXYZ)
	np.savetxt(path+'loadEnvelope/'+folderName+sectionName+'LoadEnvelopeID.txt',wlEnvID)
	np.savetxt(path+'loadEnvelope/'+folderName+sectionName+'LoadEnvelopeIDXYZ.txt',wlEnvIDXYZ)
	#makeWingLoadEnvelope('am_analysis/loads/cases/',[i for i in range(0,8)])
	
def findSparFlangesHeight(loadEnv,ix,im,spar,width,sigma,offset=0.0,k=1.0):
	x=loadEnv[:,ix]
	wingThick=np.interp(x,spar[:,0],spar[:,1])-offset
	b=BeamCrossSection()
	flangeThick=np.zeros(wingThick.shape)
	for i in range(0,flangeThick.shape[0]):
		b.setBHbh(width,wingThick[i],width,0.0)
		b.findFlangeHeightZZ(abs(loadEnv[i,im]*k),sigma)
		flangeThick[i]=(b.H-b.h)/2.
	return flangeThick

def makeSparFlangeThickness(path):
	sparThick=np.loadtxt('am_analysis/wingSparMainThickness.txt')
	sparThick[:,0]*=-1.
	sparThick[:,1]*=1e3
	sparThick=sparThick[::-1]
	l=np.loadtxt(path+'loadEnvelope/wing/wingLoadEnvelopeAbs.txt')
	t=np.zeros((l.shape[0],2))
	t[:,0]=l[:,1]
	t[:,1]=findSparFlangesHeight(l,1,6,sparThick,50,400,offset=2,k=2.25e3)
	
	np.savetxt(path+'loadEnvelope/wing/wingMainSparThickness.txt',t,header='  wing semispan y [m]	   flange thickness [mm]')
	
	
	
def examineSubsectionFlightLoad(d,WING,spar,pl,nameFolder,nameSection,incidence,dihedral,weight, side='left',ptc=None,pto=None,FLAP=None,flapOfs=0,js=0,oi=10,i1=40,i2=60):
	pts=findSubsectionSparPoints(WING,spar,j=js,oi=oi,i1=i1,i2=i2)

	fcs0=examineSubsectionForce(WING,pts)
	mms0=examineSubsectionMoment(WING,pts)
	
	cgmassW1,massW1=findSubsectionMassDistributionByAreas(WING,pl,weight)


	#if FLAP!=None:
		#flapLeftOfs=flapOfs
		#ptsFlapLeft=findSubsectionSparPoints(FLAP,spar,j=js,oi=0,i1=3,i2=10)
		#fFlapLeft=examineSubsectionForce(FLAP,ptsFlapLeft)
		#mFlapLeft=examineSubsectionMoment(FLAP,ptsFlapLeft)
		
		#for i in range(0,ptsFlapLeft.shape[0]):
			#fcs0[i+flapLeftOfs]+=fFlapLeft[i]
			#mms0[i+flapLeftOfs]+=mFlapLeft[i]

	##

	####spar_wing_l=extractPointsFromSTEPFile(d["fileWingReferenceAxis"])/1e3
	####spar_wing_r=extractPointsFromSTEPFile(d["fileWingReferenceAxis"])/1e3
	####spar_wing_r[:,1]*=-1.

	####pl=Plane3D()
	####pl.o=np.array([688.59218156,0.0,-475.835852458])*1e-3


	####pto=np.array([688.59218156,0.0,-475.835852458])*1e-3
	####WING=WING_RIGHT
	####spar=spar_wing_r
	####incidence=d['wingIncidence']
	####dihedral=5.5
	####weight=40.
	####side='right'
	####js=1;oi=10;i1=40;i2=60
	####side='reverse'


	pts=findSubsectionSparPoints(WING,spar,j=js,oi=oi,i1=i1,i2=i2)

	fcs0=examineSubsectionForce(WING,pts)
	mms0=examineSubsectionMoment(WING,pts)



	cgmassW01,massW1=findSubsectionMassDistributionByAreas(WING,pl,weight)
	cgmassW1=np.zeros((len(cgmassW01),3))
	for i in cgmassW01:
		cgmassW1[i]=cgmassW01[i]

	if side == 'reverse':
		pts=pts[::-1]
		fcs0=fcs0[::-1]
		mms0=mms0[::-1]
		cgmassW1=cgmassW1[::-1]
		massW1=massW1[::-1]

	forcesAero={}
	intForcesAero={}
	intForcesMass={}
	intForces={}
	forcesMass={}


	inm=InertiaModel()
	inm.CG=d['CG']
	inm.linearAcceleration=d['linearAcceleration']
	inm.angularAcceleration=d['angularAcceleration']
	inm.angularVelocity=d['angularVelocity']

	for i in range(0,len(fcs0)):
		forcesAero[i]=PointLoad()
		forcesAero[i].setP(pts[i])
		forcesAero[i].setF(fcs0[i])
		forcesAero[i].setM(mms0[i])
		


		forcesMass[i]=PointLoad()
		forcesMass[i].f=np.array([0,0,-100.])
		forcesMass[i].setP(cgmassW1[i])
		forcesMass[i].setForceVector(d['inertiaForce']*-1.0)
		forcesMass[i].setF(massW1[i]*-1.0*inm.getAccelerationAtPT(cgmassW1[i]))
		


		intForces[i]=PointLoad()
		intForces[i].setP(pts[i])
		intForcesAero[i]=PointLoad()
		intForcesAero[i].setP(pts[i])
		intForcesMass[i]=PointLoad()
		intForcesMass[i].setP(pts[i])



	intForcesAero[0].setM(forcesAero[0].m)
	intForcesAero[0].setF(forcesAero[0].f)
	intForcesMass[0]+=forcesMass[0]

	if ptc is not None:
		i=len(intForcesAero)
		intForces[i]=PointLoad()
		intForces[i].setP(ptc)
		intForcesAero[i]=PointLoad()
		intForcesAero[i].setP(ptc)
		intForcesMass[i]=PointLoad()
		intForcesMass[i].setP(ptc)



	if pto is not None:
		i=len(intForcesAero)
		intForces[i]=PointLoad()
		intForces[i].setP(pto)
		intForcesAero[i]=PointLoad()
		intForcesAero[i].setP(pto)
		intForcesMass[i]=PointLoad()
		intForcesMass[i].setP(pto)


	for i in range(1,len(pts)):
		intForcesAero[i]+=intForcesAero[i-1]
		intForcesAero[i]+=forcesAero[i]
		intForcesMass[i]+=intForcesMass[i-1]
		intForcesMass[i]+=forcesMass[i]

	if ptc is not None:
		i+=1
		intForcesAero[i]+=intForcesAero[i-1]
		intForcesMass[i]+=intForcesMass[i-1]


	if pto is not None:
		i+=1
		intForcesAero[i]+=intForcesAero[i-1]
		intForcesMass[i]+=intForcesMass[i-1]
		if ptc is not None:
			intForcesAero[i].m=np.copy(intForcesAero[i-1].m)
			intForcesMass[i].m=np.copy(intForcesMass[i-1].m)
				
	for i in range(0,len(intForcesAero)):
		intForces[i]+=intForcesAero[i]+intForcesMass[i]


	wingLoad=np.zeros((len(intForcesAero),9))
	wingLoadAero=np.zeros_like(wingLoad)
	wingLoadMass=np.zeros_like(wingLoad)

	wingLoadXYZ=np.zeros_like(wingLoad)
	wingLoadAeroXYZ=np.zeros_like(wingLoad)
	wingLoadMassXYZ=np.zeros_like(wingLoad)



	for i in range(0,len(intForcesAero)):
		intForces[i].csyses[0]=bc.rotY(intForces[i].csyses[0],np.zeros(3),np.deg2rad(incidence))
		intForces[i].csyses[0]=bc.rotX(intForces[i].csyses[0],np.zeros(3),np.deg2rad(dihedral))
		intForcesAero[i].csyses[0]=bc.rotY(intForcesAero[i].csyses[0],np.zeros(3),np.deg2rad(incidence))
		intForcesAero[i].csyses[0]=bc.rotX(intForcesAero[i].csyses[0],np.zeros(3),np.deg2rad(dihedral))
		intForcesMass[i].csyses[0]=bc.rotY(intForcesMass[i].csyses[0],np.zeros(3),np.deg2rad(incidence))
		intForcesMass[i].csyses[0]=bc.rotX(intForcesMass[i].csyses[0],np.zeros(3),np.deg2rad(dihedral))
		
		

		
		
		
		wingLoad[i,0:3]=intForces[i].p
		wingLoad[i,3:6]=intForces[i].localizedForce(0)
		wingLoad[i,6:9]=intForces[i].localizedMoment(0)
		
		wingLoadAero[i,0:3]=intForcesAero[i].p
		wingLoadAero[i,3:6]=intForcesAero[i].localizedForce(0)
		wingLoadAero[i,6:9]=intForcesAero[i].localizedMoment(0)
		
		wingLoadAeroXYZ[i,0:3]=intForcesAero[i].p
		wingLoadAeroXYZ[i,3:6]=intForcesAero[i].f
		wingLoadAeroXYZ[i,6:9]=intForcesAero[i].m
		
		
		
		wingLoadMass[i,0:3]=intForcesMass[i].p
		wingLoadMass[i,3:6]=intForcesMass[i].localizedForce(0)
		wingLoadMass[i,6:9]=intForcesMass[i].localizedMoment(0)
		

		
		wingLoadMassXYZ[i,0:3]=intForcesMass[i].p
		wingLoadMassXYZ[i,3:6]=intForcesMass[i].f
		wingLoadMassXYZ[i,6:9]=intForcesMass[i].f
		
		
		wingLoadXYZ[i,0:3]=intForces[i].p
		wingLoadXYZ[i,3:6]=intForces[i].f
		wingLoadXYZ[i,6:9]=intForces[i].m




##
	np.savetxt(d['path']+nameFolder+'/'+nameSection+'_LoadTotal.txt',wingLoad)
	np.savetxt(d['path']+nameFolder+'/'+nameSection+'_LoadAero.txt',wingLoadAero)
	np.savetxt(d['path']+nameFolder+'/'+nameSection+'_LoadInertia.txt',wingLoadMass)
	np.savetxt(d['path']+nameFolder+'/'+nameSection+'_LoadTotalXYZ.txt',wingLoadXYZ)
	np.savetxt(d['path']+nameFolder+'/'+nameSection+'_LoadAeroXYZ.txt',wingLoadAeroXYZ)
	np.savetxt(d['path']+nameFolder+'/'+nameSection+'_LoadInertiaXYZ.txt',wingLoadMassXYZ)

	return forcesAero, forcesMass
