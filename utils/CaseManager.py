import copy
from panelMethod.panelMethod1 import *
from panelMethod.env import GlobalEnvironment

class LoadCase(object):
	def __init__(self):
		self.pr1=VlmProblem()
		self.constants={'sp':1.0, 'rho':1.225,'chord':1.0,'span':1.0,'area':1.0}
		self.refPts=np.zeros((1,3))
		self.outputFunctions=[]
		self.geometryFile=''
		self.panelGroups={}
		self.results={}

	def _replaceListVariable(self,varDict,fList):
		for i,item in enumerate(fList):
			if type(item) == list:
				self._replaceListVariable(varDict,item)
			if type(item) == str:
				for k in varDict.keys():
					if k in item and type(varDict[k])==str and k!= item:
						fList[i]=fList[i].replace(k,varDict[k])
					elif k == item:# and type(varDict[k])!=str:
						try:
							fList[i]=eval(varDict[k])
						except:
							fList[i]=varDict[k]

	def makeGeometry(self):
		pr1=self.pr1
		with open(self.geometryFile) as f:
			code = compile(f.read(), self.geometryFile, 'exec')
			d={'pr1':pr1}
			exec(code, globals(),d)
		self.panelGroups=d['panelGroups']


	def examine(self,vd,fl,verbose=True):
		sp=self.sp
		alpha=self.alpha
		if verbose:
			print '    Processing load case : '+str(np.rad2deg(alpha))
		self.v=np.array([sp*np.cos(alpha),0.0,sp*np.sin(alpha)])
		self.pr1.setFreeVelocity(copy.copy(self.v))
		self.constants['v']=copy.copy(self.v)
		self.pr1.setRho(copy.copy(self.rho))
		self.outputFunctions=copy.deepcopy(fl)


		for i,function in enumerate(self.outputFunctions):
			ff=self.outputFunctions[i]
			self._replaceListVariable(vd,ff)
			if len(ff) == 3:
				self.results[ff[0]]=ff[1](*ff[2])
				if verbose:
					print '        processing function: '+str(ff[1])

			if len(ff) == 2:
				ff[0](*ff[1])
				if verbose:
					print '        processing function: '+str(ff[0])

class LoadCaseSequence(object):

	def __init__(self):
		self.cases={}
		self.constants={'sp':1.0, 'rho':1.225,'chord':1.0,'span':1.0,'area':1.0}
		self.results={}

	def _replaceListVariable(self,varDict,fList):
		for i,item in enumerate(fList):
			if type(item) == list:
				self._replaceListVariable(varDict,item)
			if type(item) == str:
				for k in varDict.keys():
					if k in item and type(varDict[k])==str and k!= item:
						fList[i]=fList[i].replace(k,varDict[k])
					elif k == item:# and type(varDict[k])!=str:
						try:
							fList[i]=eval(varDict[k])
						except:
							fList[i]=varDict[k]


	def makeCases(self,geometryFile,alphaRange):
		self.geometryFile=geometryFile
		self.alphaRange=alphaRange
		for i in range(0,len(alphaRange)):
			self.cases[i]=LoadCase()
			for k in self.constants:
				setattr(self.cases[i],k,self.constants[k])
			self.cases[i].alpha=alphaRange[i]
			self.cases[i].geometryFile=self.geometryFile
			self.cases[i].makeGeometry()
			self.cases[i].ID=i
			self.cases[i].constants=copy.copy(self.constants)
			self.cases[i].constants['ID']=i
			self.cases[i].constants['alphaDeg']=np.rad2deg(alphaRange[i])
			self.cases[i].constants['alpha']=alphaRange[i]
			self.cases[i].constants['geometryFile']=self.geometryFile
			for g in self.cases[i].panelGroups:
				self.cases[i].panelGroups[g].setRefPT(self.cases[i].constants['cg'])
		
	def examine(self,variables,fd,verbose=True):
		if verbose:
			print 'Processing output requests:'
		for c in self.cases:
			self.cases[c].examine(variables['caseLocal'],fd['caseLocal'])
		self.outputFunctions=copy.deepcopy(fd['global'])


		for i in range(0,len(self.outputFunctions)):
			ff=self.outputFunctions[i]
			self._replaceListVariable(variables['global'],ff)
			if len(ff) == 3:
				self.results[ff[0]]=ff[1](*ff[2])
				if verbose:
					print '    processing function: '+str(ff[1])
			if len(ff) == 2:
				ff[0](*ff[1])
				if verbose:
					print '    processing function: '+str(ff[0])

	def flushCases(self):
		del self.cases

#	def createOutputRequests(self,variables):
#		fll=copy.copy(fl)
#		replaceListVariable(variables,fll)
#		return fll
#	
#class LoadCaseManager(object):
#	def __init__(self):
#		self.constants={}
#		pass
#
#	
#	def getLoadCaseSequence(self,fl,variables):
#		lcs=LoadCaseSequence()
#		lcs.constants=copy.copy(self.constants)
#		lcs.create(copy.copy(variables))
#		return lcs
