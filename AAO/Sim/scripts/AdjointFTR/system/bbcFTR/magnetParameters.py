import sys, os, pathlib
from pathlib import Path
import numpy as np

simRoot = Path(__file__).parent.parent.parent
sys.path.append( os.path.join(simRoot,'momentSolver') )

from Magnets import QuadProfile, SolenoidProfile

###################################
# quad profile settings
g0 = 0.0361 # T / m A
l = 0.12 # meters
s0 = 0.06 # meters
d = 0.021 # meters
# solenoid profile settings
b = 5.0408
c = 0.5027
ds = 137.08
ls = 200 # cm
s0s = 100 #cm

# Quad settings
quad1 = {}
quad1['type'] = 'quad'
quad1['zstart'] = .00425
quad1['length'] = l
quad1['zend'] = quad1['zstart'] + quad1['length']
quad1['rotation'] = 45*np.pi/180
quad1['dbdx'] = QuadProfile( lambda s : -g0 * np.exp(-1 * (s-s0)**2 / d**2 ) )

quad2 = {}
quad2['type'] = 'quad'
quad2['zstart'] = 0.10655
quad2['length'] = l
quad2['zend'] = quad2['zstart'] + quad2['length']
quad2['rotation'] = 45*np.pi/180
quad2['dbdx'] = QuadProfile( lambda s : g0 * np.exp(-1 * (s-s0)**2 / d**2 ) )

quad3 = {}
quad3['type'] = 'quad'
quad3['zstart'] = 0.20895
quad3['length'] = l
quad3['zend'] = quad3['zstart'] + quad3['length']
quad3['rotation'] = 45*np.pi/180
quad3['dbdx'] = QuadProfile( lambda s : -g0 * np.exp(-1 * (s-s0)**2 / d**2 ) )

sol = {}
sol['type'] = 'solenoid'
sol['zstart'] = 0.213
sol['length'] = ls * 1e-2 # m
sol['zend'] = sol['zstart'] + sol['length']
sol['rotation'] = 0.0 # doesn't do anything
sol['dbdx'] = SolenoidProfile( lambda s,I : -1e-4*(7.4153*I - 0.12786)*c*( ((s*1e2-s0s)+ds/2)/(((s*1e2-s0s)+ds/2)**2+b**2)**0.5 - ((s*1e2-s0s)-ds/2)/(((s*1e2-s0s)-ds/2)**2+b**2)**0.5 ) )

lattice = np.array([ quad1, quad2, quad3, sol ])
###################################