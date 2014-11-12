#! /usr/bin/python
from numpy import *
import numpy
import math
import sys
from matplotlib import pyplot as plt

MAX_DENSITY=800.0
MIN_DENSITY=10.0
DEFAULT_DENSITY=500.0
### static model paramters    ###
CP = 2116.0   #specific heat capacitance of ice J/kg/K
#################################

# heat concentration
def Temp2Heat(temperature,density):
    return density*CP*temperature

def Heat2Temp(heat,density):
    return heat/(CP*density)

# heat conductivity
def K(density):
    return KSturn(density)
    #return KAbel(density)

# Abel
def KAbel(density):
    return 3.14e-6*(density**2)

# Sturn
def KSturn(density):
    return 0.138e-6+(3.233e-6*density - 1.01e-6)*density

# Exception handling classes
class IterationStagnation(Exception):
    def __init__(self,message):
        Exception.__init__(self,message)
        # NOTE: or use "super" to call base class constructor

class ExitCondition(Exception):
    def __init__(self,message):
        Exception.__init__(self,message)

# Snow node class
class SnowNode:
    # constructor
    def __init__(self, height, density, heat_concentration,inf_reservoir=array(False)):
        ### node behaviour flags  ###
        self.InfiniteReservoir = False # if True, then this node will not
                              # change in temperature, heat concentration

        ### variable model parameters ###
        self.z  = 0.0      #this snow nodes height position in the snow pack
        self.rho = 500.0   #density - kg/m3
        self.q  = None     # currently applied energy density
        self.nq = None     # new energy density
        ### linked nodes ###
        self.AboveNode = None
        self.BelowNode = None
        self.Parent = None
        ### input ###
        self.z = height
        self.rho = density
        self.q = heat_concentration
        self.InfiniteReservoir = inf_reservoir

    # print methods
    def __repr__(self):
        return "%f %f %f"%(self.z,self.T(),self.rho)

    def __str__(self):
        return "%f  %f  %f %f %s"%(self.z,self.T(),self.rho,self.q,self.InfiniteReservoir)

    # node height position in snow pack
    def Z(self):
        return self.z[0]

    def setZ(self,z):
        self.z[0] = z

    def getZ(self):
        return self.Z()

    # heat concentration
    def stepHeat(self,dq):
        if not self.InfiniteReservoir:
            self.nq[0] += dq

    # temperature
    def T(self):
        return Heat2Temp(self.q,self.rho)[0]
    
    def getT(self):
        return self.T()

    def setT(self,T):
        self.q[0] = Temp2Heat(T,self.rho)[0]

    def setCelsiusT(self,T):
        self.setT(T+273.15)

    def getCelsiusT(self):
        return self.getT()-273.15

    # density
    def getDensity(self):
        return self.rho[0]

    def setDensity(self,rho):
        if rho > 0.0:
            # preserve temperature,
            T=self.getT()
            # set density
            self.rho[0]=rho
            # reset temperature
            self.setT(T)
    
    # settings
    def setInfiniteReservoir(self,TruthValue):
        self.InfiniteReservoir[0]=TruthValue

class SnowPack:
    DEBUG_COUNTER=0

    def __init__(self,
                 heights = arange(0.0,1.0,0.05),
                 temperatures = arange(273.15,263.15,-0.5),
                 densities = None,
                 ):
        # class variables
        #################################
        self.dt = 600.0      #discrete time step - s
        self.t=0.0
        #################################
        self.SnowNodeList=[]
        self.SnowNodeDict={}
        self.TopBoundary=None
        self.BottomBoundary=None
        self.StartTemperatures=None
        ### Snow node data refs #########
        self.Heights=None
        self.Densities=None
        self.HeatConcentrations=None
        #################################

        # grid size
        N=heights.shape[0]

        # save the initial state
        self.StartTemperatures=temperatures.copy()

        if densities == None:
            densities=ones(N)*DEFAULT_DENSITY
        
        # convert temp to heat
        heatConcentrations = Temp2Heat(temperatures,densities)

        # reservoir type flags
        infiniteReservoirs = zeros(N,dtype=bool)
        
        # store numpy array reference to node data (careful to copy so not to reference external data)
        self.Heights=heights.copy()
        self.Densities=densities.copy()
        self.HeatConcentrations=heatConcentrations.copy()
        self.InfiniteReservoirs=infiniteReservoirs.copy()

        # create snow nodes
        for i in range(N):
            # create node 
            node=SnowNode(self.Heights[i:i+1],self.Densities[i:i+1],self.HeatConcentrations[i:i+1],inf_reservoir=self.InfiniteReservoirs[i:i+1])
            self.SnowNodeList.append(node)
            self.SnowNodeDict[int(round(self.Heights[i]*1000.0))]=node   # for quick lookup by height (in mm)
            # link node to previous node
            if i>0:
                self.SnowNodeList[i].BelowNode = self.SnowNodeList[i-1]
                self.SnowNodeList[i-1].AboveNode = self.SnowNodeList[i]
            # link to parent
                self.SnowNodeList[i].Parent = self

        # set boundary nodes to infinite reservoirs
        self.SnowNodeList[ 0].setInfiniteReservoir(True)
        self.SnowNodeList[-1].setInfiniteReservoir(True)


    # print routines
    def __str__(self):
        s = "--------------------------------------------------------------\n"
        s +=  "| Height | Temperature | Density | Heat conc. | Inf.Reserv. |\n"
        s += "--------------------------------------------------------------\n"
        for sn in self.SnowNodeList:
            s += str(sn) + "\n"
        s += "--------------------------------------------------------------\n"
        return s

    # configure routines
    def setDt(self,dt):
        self.dt = dt

    def getDt(self):
        return self.dt

    def setBoundary(self,topBoundary=None,bottomBoundary=None):                       
        # apply
        if topBoundary is not None:
            self.TopBoundary=topBoundary.copy()
        if bottomBoundary is not None:
            self.BottomBoundary=bottomBoundary.copy()

    # run simulation routines
    # run a single time step of simulation
    def step(self):
        # move time step forward
        self.t += self.dt
        # update boundary conditions to this step, if provided
        if self.TopBoundary is not None: 
            idx = self._find_nearest(self.TopBoundary['t'],self.t)
            self.setTopBoundTemperature(self.TopBoundary['T'][idx])

        if self.BottomBoundary is not None: 
            idx = self._find_nearest(self.BottomBoundary['t'],self.t)
            self.setBottomBoundTemperature(self.BottomBoundary['T'][idx])

        # store original heat conc. in inf reservoirs
        ir=self.InfiniteReservoirs
        infResHeat = (self.HeatConcentrations[ir]).copy()

        # calculate heat transfer
        T=self.getTemperatures()
        z=self.Heights
        dT = T[1:]-T[:-1]
        dz = z[1:]-z[:-1]
        #K  = self.getConductivities()[:-1] ## try take mean conductivity here?
        KK  = self.getConductivities()
        K = (KK[:-1]+KK[1:])/2.0
        dq = K*dT*self.dt/(dz*dz)

        # apply heat changes
        # adding changes from above
        self.HeatConcentrations[:-1] += dq
        # adding changes from below
        self.HeatConcentrations[1:] -= dq

        # concerve heat concentration in infinite reservoirs
        self.HeatConcentrations[ir] = infResHeat

    def _find_nearest(self,array,value):
        idx = (abs(array-value)).argmin()
        return idx

    # run over 'period' number of seconds 
    def run(self,stop_t):
        while self.t < stop_t:
            self.step()
            #self.DEBUG_COUNTER+=1
        #print self.DEBUG_COUNTER
        #self.DEBUG_COUNTER=0

    # return to start temperature state
    def resetToStartTemperatures(self):
        self.setTemperatures(self.StartTemperatures)
        self.t=0

    def rewind(self):
        self.resetToStartTemperatures()

    # chi-square comparison of current state with an input profile
    def chiSquare(self,ZComp,TComp):
        cs=0.0
        for i in range(len(ZComp)):
            cs += (self.getTemperatureAtHeight(ZComp[i]) - TComp[i])**2.0
        return cs/len(ZComp)

    # Entropy (as in Max. Ent.) of snow density in snow pack
    def MEMEntropy(self):
        TotRho=self.Densities.sum()
        dnorm = self.Densities/TotRho
        e = (dnorm*log(dnorm)).sum()
        return -e

    # data accessors
    def setDensities(self,densities):
        # preserver temperature / heat
        T=self.getTemperatures()
        # disallow dens<MIN_DENSITY and dens>MAX_DENSITY
        d=array(densities)
        d[(d<MIN_DENSITY)]=MIN_DENSITY
        d[(d>MAX_DENSITY)]=MAX_DENSITY
        # set new snow density
        if self.Densities is not None:
            self.Densities[:]=d
        else:
            self.Densities=d
        # set temperature / heat
        self.setTemperatures(T)

    def getDensities(self):
        return self.Densities.copy()

    def getConductivities(self):
        return K(self.Densities)

    def setHeatConcentrations(self,heatconcentrations):
        if self.HeatConcentrations is not None:
            self.HeatConcentrations[:]=array(heatconcentrations)
        else:
            self.HeatConcentrations=array(heatconcentrations)
    
    def getHeatConcentrations(self):
        return self.HeatConcentrations.copy()
 
    def setTemperatures(self,temperatures):
        self.setHeatConcentrations( Temp2Heat(temperatures,self.Densities) )

    def getTemperatures(self):
        return Heat2Temp(self.HeatConcentrations,self.Densities)

    def getCelsiusTemperatures(self):
        return self.getTemperatures()-273.15

    def getTemperatureAtHeight(self,h):
        return self.SnowNodeDict[int(h*1000.0)].getT()

    def getHeights(self):
        return self.Heights.copy()

    def setTopBoundTemperature(self,T):
        self.SnowNodeList[-1].setT(T)

    def setBottomBoundTemperature(self,T):
        self.SnowNodeList[0].setT(T)
        
class SnowPackSolver:
    # construct solver with snow pack object,
    # and target temperature profiles at one or more time steps.
    def __init__(self,snowPack,targett,targetZ,targetT):
        if abs(targett[0]-targett[1])<snowPack.dt:
            raise Exception("Exception: fit evaluation time step smaller than model time step.")
        self.SP = snowPack
        self.TargetT=targetT
        self.TargetZ=targetZ
        self.Targett=targett

    def evaluate(self,densities):
        # revert to start state
        self.SP.resetToStartTemperatures()
        # apply new density profile
        if densities is not None:
            self.SP.setDensities(densities)

        # sum up data fit over len(self.Targett) profiles
        cs = 0.0
        for i in range(1,len(self.Targett)):
            # run model forward
            self.SP.run(self.Targett[i])
            # evaluate chiSquare
            cs += self.SP.chiSquare(self.TargetZ,self.TargetT[i,:])
        cs = cs/len(self.Targett)
        # evaluate mem entropy of density profile
        s  = self.SP.MEMEntropy()
        # return result of evaluation
        return cs, s

    def solve(self, error=0.07):
        rho0 = self.SP.getDensities()
        # solve
        xopt = TopToBottomGradientSolverMEM(self.evaluate, rho0, error=error)
        print ""
        print "SNOW MODEL SOLUTION:"
        print self.SP

# custom solver
def TopToBottomGradientSolverMEM(foo, x0, error=0.07, thresh=1.0, max_step=50.0, gain=0.5, l=1.0):
    LARGE_NUMBER=9999999999.0
    MAX_STAGNATION_COUNT=5
    StagnationCount=0

    xcur = list(x0)
    L=len(x0)

    # store metric history here, for exit condition
    MHist=[]
    NCSHist=[]
    for i in range(L):
        MHist.append(-LARGE_NUMBER)
        NCSHist.append(LARGE_NUMBER)

    DropLagrange=False
    M1=-LARGE_NUMBER
    NormChiSq=LARGE_NUMBER
    NormChiSqLastLagOpt=LARGE_NUMBER

    # start iterations
    try:
        i=1
        while True:
            k = -(i%(-L))  #looks funny, but is top to bottom iteration
            
            ###
            ### Convergence
            ###
            # convergence history
            MHist.insert(0,M1)
            MHist.pop()
            NCSHist.insert(0,NormChiSq)
            NCSHist.pop()

            # reduce langrange multiplier
            if DropLagrange:
                ## reduce lagrange parameter (temporary quick fix)
                l = l/2.0
                # reset metric history
                for q in range(L):
                    MHist[q]=-LARGE_NUMBER
                # reset this trigger
                NormChiSqLastLagOpt = LARGE_NUMBER
                # reset stagnation count
                StagnationCount=0

            if NCSHist[-1]-NCSHist[0]<0.01:
                l=l*2.0
                for q in range(L):
                    NCSHist[q]=LARGE_NUMBER
                    MHist[q]=-LARGE_NUMBER
                # record this ChiSq stagnation
                StagnationCount+=1
                # exception if many consecutive stagnations (stuck in local minima)
                if StagnationCount >= MAX_STAGNATION_COUNT:
                    raise IterationStagnation("Iteration stagnated: Stuck in a local maximum.")

            if MHist[0]-MHist[-1]<0.0:
                raise IterationStagnation("Iterations stagnated before ChiSq exit condition.") 
            
            # Metric stagnation
            #if MHist[0]-MHist[-1]<0.0:
            #    raise IterationStagnation("Iterations stagnated before ChiSq exit condition.")

            ###
            ### Parameter iterations
            ###
            # point 1
            cs,s = foo(xcur)
            M1 = s - cs/(error**2.0)*l
            NormChiSq = cs/(error**2)
            
            if NormChiSqLastLagOpt == LARGE_NUMBER:
                NormChiSqLastLagOpt = NormChiSq

            # maximize MEM metric
            sys.stdout.write("\rMetric, Entropy, Chi.Sq., Lagr.Mul =%f,%f,%f,%f"%(M1,s,NormChiSq,l))
            sys.stdout.flush()
            
            # reduce lagrange multiplier condition set here
            DropLagrange = ((NormChiSqLastLagOpt-NormChiSq)>((NormChiSqLastLagOpt-1.0)/4.0))

            ## Optimize physical parameters (step 2 and 3)
            xo = xcur[k]

            # point 2
            xcur[k]+=2.0
            cs,s = foo(xcur)
            M2 = s - cs/(error**2.0)*l

            # point 3
            xcur[k]+=2.0
            cs,s = foo(xcur)
            M3 = s - cs/(error**2.0)*l

            # fit quadratic
            fac = polyfit([xo,xo+2.0,xo+4.0],[M1,M2,M3],2)
            xtarg = -0.5*fac[1]/fac[0]

            # step to improve metric (needs to be a mountain top / not valley)
            step1 = (xtarg-xo)*gain
            step2 = math.copysign(max_step,step1)

            if abs(step1)<abs(step2):
                step = step1
            else:
                step = step2
            
            if fac[0]<0.0:
                xcur[k] = xo + step
            else:
                xcur[k] = xo

            # continue iteration through physical parameters
            i += 1                

            if (cs/(error**2.0)) < thresh:
                raise ExitCondition("Exit condition reached.")

    except ExitCondition as e:
        print ""
        print e
                
    return xcur

# Analysis / data interpolation
def ResampleTempProfile(h,T,hnew):
    from scipy.interpolate import interp1d
    f=interp1d(h,T,kind='cubic',bounds_error=False)
    Tnew=f(hnew)
    Tnew[-1]=T[-1]
    return Tnew
