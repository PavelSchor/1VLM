import pandas
import numpy as np
import copy

class Derivative(object):
    def __init__(self):
        self.incr=1.

    def copy(self):
        return copy.deepcopy(self)

    def setMain(self,dom):
        self.main=dom

    def setObject(self,obj):
        self.obj=obj
    def setConstant(self,const):
        self.const=const
    def setIncrement(self,incr):
        self.incr=incr

    def getDerivatives(self):
        # print(self.obj)
        self.main.setVelocity()
        self.main.dom1.initSolution()
        self.main.fcs.inputs.initAll()
        self.main.fcs.deflectRudders()
        self.main.dom1.reinitPanelsGeometry()
        self.main.dom1.updateTimeStep()
        self.main.dom1.compute()


        f0=self.main.dom1.getForce()
        m0=self.main.dom1.getMoment()

        self.obj+=self.incr

        self.main.setVelocity()
        # print(self.main.alpha)
        # print(self.obj)
        # print(self.main.getVelocity())

        self.main.dom1.initSolution()
        self.main.fcs.deflectRudders()
        self.main.dom1.reinitPanelsGeometry()
        self.main.dom1.updateTimeStep()
        self.main.dom1.compute()

        f1=self.main.dom1.getForce()
        m1=self.main.dom1.getMoment()

        self.obj-=self.incr

        alpha=self.main.getAlpha()

        cl0=-f0[0]*np.sin(alpha)+f0[2]*np.cos(alpha)
        cd0= f0[0]*np.cos(alpha)+f0[2]*np.sin(alpha)
        cm0= m0[1]

        cl1=-f1[0]*np.sin(alpha)+f1[2]*np.cos(alpha)
        cd1= f1[0]*np.cos(alpha)+f1[2]*np.sin(alpha)
        cm1= m1[1]


        r0=np.hstack( (cd0,cl0,cm0, f0,m0))
        r1=np.hstack( (cd1,cl1,cm1, f1,m1))
        return (r1-r0)/self.incr

class TableDerivatives(object):
    def __init__(self):
        self.drv={}
        self.res={}
        self.constCol=np.ones(9,dtype=float)
        pass

    def addDerivative(self,d,name):
        self.drv[f'{name}_@_{d.incr}']=d

    def compute(self):
        for k in self.drv:
            self.res[k]=self.drv[k].getDerivatives()/self.constCol

    def getDataFrame(self):
        df=pandas.DataFrame(self.res)
        df.index=['cd','cl','cm','cfx','cfy','cfz','cmx','cmy','cmz']
        return df

