import sys, os, pathlib
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

simRoot = Path(__file__).parent.parent.parent
sys.path.append( os.path.join(simRoot,'src','AdjointFTR','momentSolver') )
sys.path.append( os.path.join(simRoot,'system','bbcFTR','modelInput') )

from MomentSolver import MomentSolver

def plotEnvelope(momObj):

    xr = momObj.y[0,:] + momObj.y[1,:]
    yr = momObj.y[0,:] - momObj.y[1,:]
    plt.figure()
    plt.plot(momObj.z,xr, label='xrms')
    plt.plot(momObj.z,yr, label='yrms')
    plt.plot(momObj.z,momObj.ksol * 1e-6)
    plt.plot(momObj.z,momObj.kquad * 1e-8)
    plt.grid(True)
    plt.xlabel('Z position [m]')
    plt.ylabel('Moments [m]')
    plt.legend()

def setRestrictions(momObj, paramArray):

    if True:
        for i,param in enumerate(paramArray[0:-1]):
            paramArray[i] = np.clip(param, 0, 10)

    if True:
        #length = 0.12 # 12 cm
        #quad3Center = momObj.lattice[2]['zstart'] * paramArray[6] + momObj.lattice[2]['length'] / 2.0
        #paramArray[9] = (length + quad3Center) / (momObj.latticeDefault[-1]['zstart'])
        paramArray[9] = 0.8
        #print( "POS: " + str(quad3Center) + "|" + str(paramArray[9] * momObj.lattice[3]['zstart']) )

    return paramArray

def plotParams(params):
    aa = np.zeros(( len(params),len(params[0]) ))
    for i,param in enumerate(params):
        aa[i,:] = param
    
    plt.figure()
    plt.plot(aa)

###################################
# Magnet parameters
from magnetParameters import lattice
###################################

###################################
# Opt parameters
from optParameters import params,getParamArray,getParamObj
###################################

# physics settings
energy = 5e3 # eV
current = 1.0e-3 # Amps
pipeRadius = 0.0 # meters

# sim parameters
zInterval = (0, 0.7)
stepSize = 0.0001

print('Setting up initial lattice')
mom = MomentSolver(lattice, energy=energy, current=current, pipeRadius=pipeRadius, zInterval=zInterval, stepSize=stepSize)
mom.UpdateLattice(params = params)

# run moments and adjoint equations
print('Running Mom. Eqn.')
mom.Run()
mom.RunAdjoint()

# plot before
plotEnvelope(mom)

# get FoM
print('Starting Opt.')
f0,f0p,_ = mom.GetFoM_And_DFoM()
df0 = mom.GetDF()
gamma = f0 / np.sum( df0**2 )

# opt history
an_h = [getParamArray(params)]
gamma_h = [gamma]
f_h = [f0]
fp_h = [f0p]
df_h = [df0]

# initial first step
an_h.append( an_h[0] - gamma_h[0] * df_h[0] )
mom.UpdateLattice( params = getParamObj(an_h[-1]) )
an_h[-1] = setRestrictions(mom, an_h[-1])
mom.UpdateLattice( params = getParamObj(an_h[-1]) )
mom.Run()
ftmp,fptmp,_ = mom.GetFoM_And_DFoM()
f_h.append(ftmp)
fp_h.append(fptmp)
print('FoM: ' + str(f_h[-1]))

# find the starting gamma value
while f_h[-1] >= f0:
    gamma_h.append( gamma_h[-1] / 2.0 )
    an_h.append( an_h[0] - gamma_h[-1] * df_h[0] )
    an_h[-1] = setRestrictions(mom, an_h[-1])    
    mom.UpdateLattice( params = getParamObj(an_h[-1]) )
    mom.Run()
    ftmp,fptmp,_ = mom.GetFoM_And_DFoM()
    f_h.append(ftmp)
    fp_h.append(fptmp)
    print('FoM: ' + str(f_h[-1]))

# main loop 
try:
    while True:

        # step
        ii=1
        while f_h[-1] < f_h[-2]:
            print('Iterating ' + str(ii))

            # iterate
            an_h.append( an_h[-1] - gamma_h[-1] * df_h[-1] )       
            an_h[-1] = setRestrictions(mom, an_h[-1])        
            mom.UpdateLattice( params = getParamObj(an_h[-1]) )
            mom.Run()
            ftmp,fptmp,_ = mom.GetFoM_And_DFoM()
            f_h.append(ftmp)
            fp_h.append(fptmp)
            print('FoM: ' + str(f_h[-1]))
            ii += 1

            # if we have done 20 steps in a row, start increasing the step size
            if( ii > 20 ):
                gamma_h.append( gamma_h[-1] * 2.0 )

        # can't step anymore, recompute adjoint
        print('Recomputing Adjoint Equations')

        # grab last good setting
        an_h.append( an_h[-2] )
        f_h.append( f_h[-2] )
        fp_h.append( fp_h[-2] )

        # calculate adjoint
        mom.UpdateLattice( params = getParamObj(an_h[-1]) )
        mom.Run()
        mom.RunAdjoint()

        # calculate df
        df_h.append( mom.GetDF() )

        # no improvement from last step, try to update gamma
        if( ii == 2 ):
            print('Updating Gamma')
            f0n = f_h[-1]
            ann = an_h[-1]

            iii = 1
            while f_h[-1] >= f0n:
                gamma_h.append( gamma_h[-1] / 2.0 )           
                an_h.append( ann - gamma_h[-1] * df_h[-1] )
                an_h[-1] = setRestrictions(mom, an_h[-1])            
                mom.UpdateLattice( params = getParamObj(an_h[-1]) )
                mom.Run()
                ftmp,fptmp,_ = mom.GetFoM_And_DFoM()
                f_h.append(ftmp)
                fp_h.append(fptmp)
                print('FoM: ' + str(f_h[-1]))
                iii += 1

                if ( iii > 25):
                    break

        if ( f_h[-1] < 1e-14 ):
            break

        if ( len(f_h) > 1000 ):
            break

        if ( ii == 1 ):
            break
except KeyboardInterrupt:
    pass

print(an_h[-1])
plotEnvelope(mom)

plt.show()