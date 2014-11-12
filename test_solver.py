#! /usr/bin/python
from sys import exit, argv
from os.path import basename
import matplotlib.pyplot as plt
from numpy import array, arange, insert, ones, around

# Read data
def ReadData(fnames):
    TEMPSTEP=0.03125
    times=[]
    Temps=[]
    for fname in fnames:
        f=open(fname,"r")
        Lines=f.readlines()
        for l in Lines:
            splt=l.split(':') 
            times.append(float(splt.pop(0)))
            Tprof=[]
            for T in splt:
                Tprof.insert(0,float(T)*TEMPSTEP)
            Temps.append(Tprof)
        f.close()
    return array(times), arange(0.2,3.4,0.2), array(Temps)+273.15

# Read data
argv.pop(0)
t,h,T=ReadData(argv)
t=t-t[0]
t=around(t/600.0)*600.0  #round to whole 600.0s , because snow pack does not yet handle irregular time spacing.

# chose time offset
tidx=len(t)-1

# select initial state data
sel=(h<1.81)
hprof=h[sel]
Tprof=(T[0,:])[sel]

# Snow Pack simulation
from simpssim import SnowPack, SnowPackSolver, ResampleTempProfile, IterationStagnation
from numpy import random

##### SETTINGS ##################
hinit=arange(0.2,1.81,0.2)
Tinit=ResampleTempProfile(hprof,Tprof,hinit)

N_SOLVES = 4 # number of independent solutions (good to check variance in solution)
SOLVER_TIME_STEP = 1200.0 # use multiples of SM4 time step 600s
DATA_ERROR = 0.07  # uniform error in temperature data
################################

Dinitial=[]
Dresults=[]
for i in range(N_SOLVES):
    # add perturbation to initial temperature condition
    error = DATA_ERROR
    pertTinit = Tinit# + random.normal(0.0,error,len(Tinit))
    # random initial condition
    dinit=hinit[::-1]*0.0+100.0+random.rand(len(hinit))*599.0
    sp=SnowPack(heights=hinit,temperatures=pertTinit,densities=dinit)
    sp.setDt(SOLVER_TIME_STEP)
    print "INITIAL SNOW MODEL CONDITIONS:"
    print sp

    # boundary condition
    sp.setBoundary(topBoundary={'t':t[:],'T':T[:,8]},bottomBoundary={'t':t[:],'T':T[:,0]})

    # set up solver
    lt = t.shape[0]

    spsL=SnowPackSolver(sp,t[::6],h[sel],T[::6,sel])
    try:
        print "Solving: %d/%d"%(i+1,N_SOLVES)
        spsL.solve(error=error)
    except IterationStagnation as e:
        print ""
        print e
    else:
        Dresults.append(sp.getDensities())
        Dinitial.append(dinit)
    
# show variation in initial profiles
datestr = basename(argv[0])[4:6]+"/"+basename(argv[0])[2:4]+"/20"+basename(argv[0])[0:2]
for d in Dinitial:
    plt.plot(d,sp.getHeights(),'r--')
    plt.plot(d,sp.getHeights(),'bx')
plt.title("snow density start conditions\n before optimization")
plt.xlabel("snow density (kg/m3)")
plt.ylabel("height (m)")
plt.xlim((0.0,800.0))
plt.ylim((0.0,2.0))
plt.savefig('density_initial.png')
plt.show()


# show inferred density profiles
for d in Dresults:
    plt.plot(d,sp.getHeights(),'r--')
    plt.plot(d,sp.getHeights(),'bx')
plt.title("snow density solutions to SM4 profiles\n covering period %s 00:00 + %d hours"%(datestr,len(argv)*24))
plt.xlabel("snow density (kg/m3)")
plt.ylabel("height (m)")
plt.xlim((0.0,800.0))
plt.ylim((0.0,2.0))
plt.savefig('density_solutions.png')
plt.show()


# calculate stdev of mean densities
Dav=0.0
Dst=0.0
for d in Dresults:
    Dav += d.sum()/len(d)
Dav=Dav/len(Dresults)
for d in Dresults:
    Dst += (Dav - d.sum()/len(d))**2.0
Dst = (Dst/len(Dresults))**0.5

print "Average 'mean density' of solutions + stdev",Dav,Dst

# show animation
sp.resetToStartTemperatures()
sp.setDt(600.0)
plt.ion()
fig = plt.figure()
ax  = fig.add_subplot(111)
ax.set_xlim((-6,0.0))
ax.set_ylim((0.0,3.4))
ax.set_xlabel("temperature (deg.C)")
ax.set_ylabel("height (m)")
lineL, = ax.plot(sp.getCelsiusTemperatures(),sp.getHeights(), 'b-', label="simulation/fit")
bmarks, = ax.plot(sp.getCelsiusTemperatures()[[0,-1]],sp.getHeights()[[0,-1]], 'bo', label="boundary cond.")
xmarks, = ax.plot(T[0,:]-273.15, h, 'rx', markersize=10, label="measurements")
ax.legend()
import time

for i in range(0,len(t),1):
    sp.run(t[i])
    lineL.set_data(sp.getCelsiusTemperatures(),sp.getHeights())
    bmarks.set_data(sp.getCelsiusTemperatures()[[0,-1]],sp.getHeights()[[0,-1]])
    xmarks.set_data(T[i,:]-273.15,h)
    ax.set_title("snow pack solution compared to SM4 profiles\n time step %s + %02d:%02d"
                 %(datestr,t[i]/3600.0,(t[i]/60.0)%60.0))
    fig.canvas.draw()
    plt.draw()
    fig.savefig("temp_profile_fit_%03d.png"%(i))
    time.sleep(0.25)

plt.show()




