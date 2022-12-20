import sys, os, pathlib
from pathlib import Path
import numpy as np

#simRoot = Path(__file__).parent.parent.parent
#sys.path.append( os.path.join(simRoot,'src','AdjointFTR','momentSolver') )
#sys.path.append( os.path.join(simRoot,'system','bbcFTR','modelInput') )

from momentSolver.MomentSolver import MomentSolver
import numpy as np
import matplotlib.pyplot as plt

###################################
# Magnet parameters
from system.bbcFTR.magnetParameters import lattice
###################################

###################################
# Opt parameters
from system.bbcFTR.optParameters import params
###################################

# physics settings
energy = 5e3 # eV
current = 0.0e-3 # Amps
pipeRadius = 0.0 # meters

# sim parameters
zInterval = (0, 0.322)
stepSize = 0.0001

mom = MomentSolver(lattice, energy=energy, current=current, pipeRadius=pipeRadius, zInterval=zInterval, stepSize=stepSize)
mom.optParams = params
#mom.UpdateLattice(params = params)

z, y, ksol, kquad = mom.Run()
zadj, yadj, _, _ = mom.RunAdjoint()

xr = y[0,:] + y[1,:]
yr = y[0,:] - y[1,:]
plt.figure()
plt.plot(z,xr)
plt.plot(z,yr)
plt.plot(z,ksol * 1e-7)
plt.plot(z,kquad * 1e-9)
plt.grid(True)

plt.show()