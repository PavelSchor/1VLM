import numpy as np
import copy

from utils.b3Vect import B3vect
import utils.bCad

from vlm.vlmPanels import VlmPanel

bc=utils.bCad.BCad()
vect=B3vect()

def getReversedPanelGroup(pg,rk):
	ng=PanelGroup()
	kk=pg.panels.keys()[::-1]
	for i in kk:
		p=pg.getPanel(i)
		###!!rk=p.pts[2],p.pts[3],p.pts[0],p.pts[1])
		n=doubletSourcePanel(p.pts[rk[0]],p.pts[rk[1]],p.pts[rk[2]],p.pts[rk[3]])
		n.isWake=False
		n.connectivity=copy.deepcopy(p.connectivity)
		n.connectivity['d1']=p.connectivity['d3']
		n.connectivity['d3']=p.connectivity['d1']
		n.gid=p.gid
		if p.isFirst:
			n.isLast=True
		if p.isLast:
			n.isFirst=True
		ng.addPanel(n)
	return ng

def extractGidsByConnectivity(pr1,po,cn,checkEnds=True):
	g=[po]
	p=pr1.dom1.getPanelByGid(po)
	if checkEnds:
		while not ( p.connectivity[cn] is None or p.isLast or p.isFirst):
			p=pr1.dom1.getPanelByGid(p.connectivity[cn])
			#if not p.connectivity[cn] is None:
			g.append(p.gid)
	else:
		while not ( p.connectivity[cn] is None):
			p=pr1.dom1.getPanelByGid(p.connectivity[cn])
			#print p.gid
			#if not p.connectivity[cn] is None:
			g.append(p.gid)
	try:
		p=pr1.dom1.getPanelByGid(p.connectivity[cn])
		g.append(p.gid)
	except:
		pass
	return g

def extractLineByGids(pr1,g,j):
	pts=np.zeros((len(g),3))
	for i in range(0,len(g)):
		p=pr1.dom1.getPanelByGid(g[i])
		pts[i]=p.pts[j]
	return pts



def commonEdgePointsTwoPanels(A,B):
	nrows, ncols = A.pts.shape
	dtype={'names':['f{}'.format(i) for i in range(ncols)], 'formats':ncols * [A.pts.dtype]}
	C = np.intersect1d(A.pts.view(dtype), B.pts.view(dtype))
	if commonEdgeTwoPanels(A,B):
		r=np.zeros((2,3))
		for i in range(0,3):
			r[0,i]=C[0][i]
			r[1,i]=C[1][i]
		return r
		#return C
	else:
		return None


def commonEdgeTwoPanels(A,B):
	nrows, ncols = A.pts.shape
	dtype={'names':['f{}'.format(i) for i in range(ncols)], 'formats':ncols * [A.pts.dtype]}
	C = np.intersect1d(A.pts.view(dtype), B.pts.view(dtype))
	if len(C) < 2:
		return False
	else:
		return True



def gridLoftedPlate(spp,chs,p1,p2,p3,p4):
	v1=p4-p1
	v2=p3-p2
	v12=p2-p1
	v43=p4-p3

	spanSpacing=np.linspace(0,1,spp+1)
	chordSpacing=np.linspace(0,1,chs+1)
	grid=[]
	
	VLE=p2-p1
	VTE=p3-p4
	for i in range(0,spp):
		grid.append([])
		LE1=p1+VLE*spanSpacing[i]
		LE2=p1+VLE*spanSpacing[i+1]
		
		TE1=p4+VTE*spanSpacing[i]
		TE2=p4+VTE*spanSpacing[i+1]
		
		#print i,spanSpacing[i],LE1,LE2,TE1,TE2
		for j in range(0,chs):
			pp1=LE1+chordSpacing[j]*(TE1-LE1)
			pp2=LE2+chordSpacing[j]*(TE2-LE2)
			pp3=LE2+chordSpacing[j+1]*(TE2-LE2)
			pp4=LE1+chordSpacing[j+1]*(TE1-LE1)
			grid[i].append([])
			grid[i][j]=[pp1,pp2,pp3,pp4]
			
	return grid


def makeSpanwiseLoftedPlate(pr1,spp,chs,p1,p2,p3,p4):

	spanWiseGrid=gridLoftedPlate(spp,2,p1,p2,p3,p4)
	
	for i in range(0,len(spanWiseGrid[0])):
		region=subSection()
		ptsS=spanWiseGrid[0][i]#[0]
		chordWiseGrid=gridLoftedPlate(2,chs,ptsS[0],ptsS[1],ptsS[2],ptsS[3])
		v14=vect.makeVect(ptsS[0],ptsS[3])
		v23=vect.makeVect(ptsS[1],ptsS[2])
		region.refPT1=ptsS[0]+0.25*v14
		region.refPT2=ptsS[1]+0.25*v23

		for j in range(0,len(chordWiseGrid)):
			ptsC=chordWiseGrid[j][0]#[j]
			panel=VlmPanel(ptsC[0],ptsC[1],ptsC[2],ptsC[3])
			panel.makeVortexLines()
			panel.makeWake()
			region.addPanel(panel)
			
		pr1.addRegion(region)

def findNodeOrder(pts,ix=0,iy=1):
	
	idx=np.array([0,1,2,3],dtype=int)
	le=pts[pts[:,ix].argsort()][0:2]
	lei=idx[pts[:,ix].argsort()][0:2]

	i1,i2=lei[le[:,iy].argsort()]
	ii3,ii4=np.delete(idx,[i1,i2])

	mpt=pts.sum(axis=0)/4.
	vle=pts[i2]-pts[i1]
	v1=pts[i1]-mpt
	n1=np.cross(vle,v1)

	vte=pts[ii4]-pts[ii3]
	v3=pts[ii3]-mpt
	n3=np.cross(vte,v3)
	
	if np.dot(n1,n3) >0:
		i3=ii3
		i4=ii4
	else:
		i3=ii4
		i4=ii3
		
	
	return [i1,i2,i3,i4]



def findNodeOrder2(pts,ix=0,iy=1):
	##le=lastEdge
	idx=np.array([0,1,2,3],dtype=int)
	le=pts[pts[:,iy].argsort()][0:2]
	lei=idx[pts[:,iy].argsort()][0:2]

	i1,i4=lei[le[:,ix].argsort()]
	ii2,ii3=np.delete(idx,[i1,i4])

	mpt=pts.sum(axis=0)/4.
	vbe=pts[i1]-pts[i4]
	v1=pts[i1]-mpt
	n1=np.cross(vbe,v1)

	vte=pts[ii3]-pts[ii2]
	v3=pts[ii3]-mpt
	n3=np.cross(vte,v3)
	
	if np.dot(n1,n3) >0:
		i2=ii2
		i3=ii3
	else:
		i2=ii3
		i3=ii2
		
	
	return [i1,i2,i3,i4]




def readSalomeDAT(fname):

	with open(fname,'r') as f:
		nNodes,nElem=[int(x) for x in f.readline().split()]

	#nodes=np.loadtxt(fname,skiprows=1,max_rows=nNodes,usecols=(1,2,3))
	nd={}
	ed={}
	et={}
	quadsIds={}
	nodalConn={}
	elConn={}
	elConn0={}
	midPts={}
	with open(fname,'r') as f:
		ll=f.readlines()
		for i in range(1,nNodes+1):
			l=ll[i].split()#[int(x) for x in  ll[i].split()]
			nd[int(l[0])]=np.array([ float(l[1]),float(l[2]),float(l[3])])
		
			
		for i in range(1+nNodes,nNodes+nElem+1):
			l=[int(x) for x in  ll[i].split()][0]
			nodalConn[l]=[]
		for i in range(1+nNodes,nNodes+nElem+1):
			l=[int(x) for x in  ll[i].split()]
			eid,et[eid],nids=l[0],l[1],l[2::]
			ed[eid]=nids
			for j in nids:
				nodalConn[j].append(eid)
			elConn[eid]=[]
			elConn0[eid]=[]
		
		for eid in ed:
			for nid in ed[eid]:
				ndConn=nodalConn[nid]
				for k in ndConn:
					if et[k]==204:
						elConn0[eid].append(k)
		for eid in elConn0:
			if et[eid]==204:
				values=np.array(elConn0[eid],dtype=int)
				for k in elConn0[eid]:
					ii = np.where(values == k)[0]
					if len(ii)>1:
						elConn[eid].append(k)
				s=set(elConn[eid])
				s.remove(eid)
				elConn[eid]=list(s)
			
				midPts[eid]=np.array([nd[j] for j in ed[eid] ]).sum(axis=0)/len(ed[eid])
	return nd,ed,et,elConn,nodalConn,midPts

def buildElementOrder(fname,firstPanel,xDir,yDir,nChordwise,nSpanwise):
	nd,ed,et,elConn,nodalConn,midPts=readSalomeDAT(fname)
	#xDir=np.array([1,0,0]);yDir=np.array([0,0,1])
	#nChordwise=22
	#nSpanwise=10
	#firstPanel=317
	dpX={}
	dpY={}
	for eid in elConn:
		dpX[eid]=np.zeros(len(elConn[eid]))
		dpY[eid]=np.zeros(len(elConn[eid]))
		for k in range(0,len(elConn[eid])):
			v=midPts[elConn[eid][k]]-midPts[eid]
			v/=np.linalg.norm(v)
			dpX[eid][k]=np.dot(v,xDir)
			dpY[eid][k]=np.dot(v,yDir)
	elOrder=np.zeros(len(midPts),dtype=np.int)
	leadingEdge=firstPanel
	eid=firstPanel
	elOrder[0]=firstPanel
	kk=1
	for s in range(0,nSpanwise-1):
		for c in range(0,nChordwise-1):
			eid=elConn[eid][dpX[eid].argmax()]
			elOrder[kk]=eid
			kk+=1
		leadingEdge=elConn[leadingEdge][dpY[leadingEdge].argmax()]
		eid=leadingEdge
		elOrder[kk]=eid
		kk+=1
	for c in range(0,nChordwise-1):
		eid=elConn[eid][dpX[eid].argmax()]
		elOrder[kk]=eid
		kk+=1
	return elOrder

def extractQuadsSalomeDAT(fname,elOrder,ix=0,iy=1):
	nd,ed,et,elConn,nodalConn,midPts=readSalomeDAT(fname)
	quadsIds={}
	for i in et:
		if et[i]==204:
				quadsIds[i]=ed[i]

	quads={}
	for i in quadsIds:
		quads[i]=np.array([nd[j] for j in ed[i] ])#[nodeOrder]

	quadsVLM=np.zeros((len(ed)+1,4,3))
	for i in quads:
		nodeOrder=findNodeOrder2(quads[i],ix,iy)
		#print(nodeOrder)
		quadsVLM[i]=quads[i][nodeOrder]

	return quadsVLM[elOrder]

def makeVLMRegionSalomeDAT(pr1,region,fname,elOrder,nChordwise,nSpanwise,ix=0,iy=1,scl=1):
	quads=extractQuadsSalomeDAT(fname,elOrder,ix,iy)*scl
	
	pcnt=0
	for i in range(0,len(quads)):
		panel=VlmPanel(quads[i][0],quads[i][1],quads[i][2],quads[i][3])
		panel.pgid=pcnt
		if (i+1)%(nChordwise)==0:
			panel.isLast=True
			
		try:
			panel.gid=pr1.dom1.panelCounter
			pr1.dom1.panelCounterIncrease()
		except:
			pass
		
		pcnt+=1
		region.addPanel(panel)
	region.nSpanwise=nSpanwise
	region.nChordwise=nChordwise
	return region

def makeSpanwiseLoftedPlate1Reg(pr1,region,spp,chs,p1,p2,p3,p4):

	#spanWiseGrid=gridLoftedPlate(spp,2,p1,p2,p3,p4)
	grid=gridLoftedPlate(spp,chs,p1,p2,p3,p4)
	pcnt=0
	region.nSpanwise=spp
	region.nChordwise=chs
	for i in range(0,len(grid)):
		for j in range(0,len(grid[0])):
			ptsC=grid[i][j]
			panel=VlmPanel(ptsC[0],ptsC[1],ptsC[2],ptsC[3])
			panel.pgid=pcnt
			
			try:
				panel.gid=pr1.dom1.panelCounter
				pr1.dom1.panelCounterIncrease()
			except:
				pass
			
			pcnt+=1
			region.addPanel(panel)

def makeMeshFromArrayList(pr1,region,ll,nip=None):
	if nip is None: nip=np.ones(len(ll)-1,dtype=int)
	if isinstance(nip,list): nip=np.array(nip,dtype=int)
	ns=nip.sum()
	nc=len(ll[0])-1
	pctn=0
	region.nSpanwise=ns
	region.nChordwise=nc
	
	
	for i in range(1,len(ll)):
		g1=ll[i-1]
		g2=ll[i]
		v=g2-g1
		rng=np.linspace(0,1,int(nip[i-1]+1))
		for j in range(0,len(rng)-1):
			s1=g1+v*rng[j]
			s2=g1+v*rng[j+1]
			for k in range(0,nc):
				panel=VlmPanel(s1[k],s2[k],s2[k+1],s1[k+1])
				panel.pgid=pctn
				if k == nc-1: panel.isLast=True
				try:
					panel.gid=pr1.dom1.panelCounter
					pr1.dom1.panelCounterIncrease()
				except:
					pass
				
				pctn+=1
				region.addPanel(panel)
					

def makeSpanwiseLofts(pr1,region,sts):
	m=sts.shape[0]
	n=sts.shape[1]
	pcnt=0
	region.nSpanwise=m
	region.nChordwise=n
	for i in range(1,m):
		for j in range(1,n):
			panel=VlmPanel(sts[i-1][j-1],sts[i][j-1],sts[i][j],sts[i-1][j])
			panel.pgid=pcnt
			
			try:
				panel.gid=pr1.dom1.panelCounter
				pr1.dom1.panelCounterIncrease()
			except:
				pass
			
			pcnt+=1
			region.addPanel(panel)



def appendLoftedWakeRegKnownLasts(pr1,region,chs,wakeVector,i2=2,i3=3):###,p1,p2,p3,p4):
	spp=region.nSpanwise
	chsr=region.nChordwise
	for pId in region.panels:
		pp0=region.panels[pId]
		if pp0.isLast:
			lastPanel=pp0
			p1=pp0.pts[i3]
			p2=pp0.pts[i2]
			p3=p2+wakeVector
			p4=p1+wakeVector
			grid=gridLoftedPlate(1,chs,p1,p2,p3,p4)
			pcnt=0
			spp=region.nSpanwise
			
			
			lastPanel.isLast=True
			lastPanel.wakes={}
			gg=lastPanel.gid
			#print(gg,lastPanel)
			region.wakes[gg]={}
			region.wakes[gg]['panels']={}
			for j in range(0,chs):
				ptsC=grid[0][j]
				region.wakes[gg]['panels'][j]=VlmPanel(ptsC[0],ptsC[1],ptsC[2],ptsC[3])
				#panel.pgid=pcnt
				lastPanel.wakes[j]=region.wakes[gg]['panels'][j]
				#try:
					#panel.gid=pr1.dom1.panelCounter
					#pr1.dom1.panelCounterIncrease()
				#except:
					#pass
				print(lastPanel.wakes)
			#pcnt+=1
			#region.addPanel(panel)

def appendLoftedWakeReg(pr1,region,chs,wakeVector):###,p1,p2,p3,p4):
	spp=region.nSpanwise
	chsr=region.nChordwise
	
	pp0=region.getPanelByPGID(chsr-1)
	pp1=region.getPanelByPGID(region.nPanels-1)
	p1=pp0.pts[3]
	p2=pp1.pts[2]
	p3=p2+wakeVector
	p4=p1+wakeVector
	#spanWiseGrid=gridLoftedPlate(spp,2,p1,p2,p3,p4)
	grid=gridLoftedPlate(spp,chs,p1,p2,p3,p4)
	pcnt=0
	spp=region.nSpanwise
	
	
	for i in range(0,spp):
		lastPanel=region.getPanelByPGID( chsr-1 +i*chsr )
		lastPanel.isLast=True
		lastPanel.wakes={}
		gg=lastPanel.gid
		print(gg,lastPanel)
		region.wakes[gg]={}
		region.wakes[gg]={}
		region.wakes[gg]['panels']={}
		for j in range(0,chs):
			ptsC=grid[i][j]
			region.wakes[gg]['panels'][j]=VlmPanel(ptsC[0],ptsC[1],ptsC[2],ptsC[3])
			#panel.pgid=pcnt
			lastPanel.wakes[j]=region.wakes[gg]['panels'][j]
			#try:
				#panel.gid=pr1.dom1.panelCounter
				#pr1.dom1.panelCounterIncrease()
			#except:
				#pass
			print(lastPanel.wakes)
			#pcnt+=1
			#region.addPanel(panel)
