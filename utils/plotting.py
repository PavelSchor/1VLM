import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
sns.set_style("whitegrid")
import json_tricks as json
import matplotlib.lines as lines



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
		self.showLegend=True
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
		
		#plt.rc('legend',fontsize=30)
		fig = plt.figure()
		ax = plt.subplot(111)


		for l in self.lines:
			
			ax.plot(l.x,l.y,label=l.label,linewidth=l.linewidth,marker = self.markers[l.ID%len(self.markers)],markersize=5,markevery=5)


		plt.grid(True)
		plt.xlabel(self.xlabel)
		plt.ylabel(self.ylabel)
		plt.title(self.title)
		if self.showLegend:
			# Shrink current axis by 20%
			#box = ax.get_position()
			#ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])

			## Put a legend to the right of the current axis
			#ax.legend(loc='upper left', bbox_to_anchor=(0.95, 0.9),ncol=1,mode='expand')
			leg=ax.legend(loc='upper right',ncol=1,framealpha=0.)#,mode='expand')
			plt.setp(plt.gca().get_legend().get_texts(), fontsize='12')
			leg.get_frame().set_alpha(1.)
		if not fname is None:
			plt.savefig(fname, dpi=400)
		else:
			plt.show()
		plt.close()	

