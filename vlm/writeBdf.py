import collections
from fem.pyNastran.bdf.field_writer_8 import print_card_8


def unique(a):
	order = np.lexsort(a.T)
	a = a[order]
	diff = np.diff(a, axis=0)
	ui = np.ones(len(a), 'bool')
	ui[1:] = (diff != 0).any(axis=1) 
	return a[ui]


gidHT={}
conn={}
gids=pr1.dom1.gids.keys()
j=0
for i in range(0,len(gids)):
	p=pr1.dom1.getPanelByGid(gids[i])
	if not p.isWake:
		gidHT[j]=gids[i]
		conn[j]=np.array([0,1,2,3],dtype=np.int)+(4*i)
		j+=1
		
		
grid0=np.zeros(( len(gidHT)*4,3))

for i in gidHT:
	p=pr1.dom1.getPanelByGid(gidHT[i])
	p.eid=i+1
	grid0[i*4:i*4+4,:]=p.pts

#grid=unique(grid0)
grid=grid0

elem=collections.OrderedDict()

eid=1
for i in gidHT:
	p=pr1.dom1.getPanelByGid(gidHT[i])
	elem[eid]={}
	elem[eid]['EID']=eid
	elem[eid]['PID']=1
	elem[eid]['DP']=p.getPressureCoef()*0.5*p.getRho()*np.linalg.norm(pr1.getFreeVelocity())**2.
	etype='CQUAD4'

	enodes=unique(p.pts)
	if len(enodes) == 3:
		etype='CTRIA3'

	elem[eid]['etype']=etype

	if etype=='CQUAD4':
		elem[eid]['conn']=[0,0,0,0]
		for j in range(0,4):
			elem[eid]['conn'][j]=1+np.where(np.all(p.pts[j] == grid,axis=1))[0][0]


	if etype=='CTRIA3':
		elem[eid]['conn']=[0,0,0]
		for j in range(0,3):
			elem[eid]['conn'][j]=1+np.where(np.all(p.pts[j] == grid,axis=1))[0][0]
	eid+=1

f=open('000pm.bdf','w')
sc=1e3
f.write('''SOL 101
CEND
ECHO = NONE
SUBCASE 1
   SUBTITLE=Default
   SPC = 2
   LOAD = 2
   DISPLACEMENT(SORT1,REAL)=ALL
   SPCFORCES(SORT1,REAL)=ALL
   OLOAD(SORT1,REAL)=ALL
   ESE=ALL
   STRAIN(SORT1,REAL,VONMISES,STRCUR,BILIN)=ALL
   GPFORCE=ALL
   STRESS(SORT1,REAL,VONMISES,BILIN)=ALL
   FORCE(SORT1,REAL,BILIN)=ALL
$ Direct Text Input for this Subcase
BEGIN BULK
MDLPRM,HDF5,1
PARAM    POST    0
PARAM   PRTMAXIM YES''')
f.write('\n')



for e in elem:
	eid=elem[e]['EID']
	pid=elem[e]['PID']
	etype=elem[e]['etype']
	cn=elem[e]['conn']
	if etype=='CQUAD4':
		f.write('{:8s}{:8d}{:8d}{:8d}{:8d}{:8d}{:8d}\n'.format('CQUAD4', eid, pid ,cn[0],cn[1],cn[2],cn[3]))

	if etype=='CTRIA3':
		f.write('{:8s}{:8d}{:8d}{:8d}{:8d}{:8d}\n'.format('CTRIA3', eid, pid ,cn[0],cn[1],cn[2]))


for i in range(0,grid.shape[0]):
	f.write('{:8s}{:16d}{:16s}{:16.8f}{:16.8f}{:8s}\n{:8s}{:16.8f}{:16s}{:16s}{:16s}{:8s}\n'.format('GRID*', i+1, '',sc*grid[i,0],sc*grid[i,1],'','*',sc*grid[i,2],'','','',''))

f.write('LOAD     2      1.      1.       1\n')
for e in elem:
	eid=elem[e]['EID']
	dp=elem[e]['DP']
	#if etype=='CQUAD4':
		#f.write('{:8s}{:8d}{:16.8f}{:8d}\n'.format('PLOAD2', 2, dp ,eid))
	f.write(print_card_8(['PLOAD2', 2,dp,e]))
f.close()