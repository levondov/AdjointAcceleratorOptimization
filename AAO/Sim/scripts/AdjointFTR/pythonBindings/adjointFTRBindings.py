import ctypes
import os
import sys, os, pathlib
from pathlib import Path
import numpy as np

simRoot = Path(__file__).parent.parent.parent.parent.parent.parent
dllRoot = os.path.join(simRoot,'out','build','x64-release','AAO','Sim','src','AdjointFTR')

lib = ctypes.WinDLL(dllRoot + os.sep + 'AdjointFTR.dll')

lib.AdjointFTR_new.restype = ctypes.POINTER(ctypes.c_char)
lib.AdjointFTR_getSCVM.argtypes = (
	ctypes.POINTER(ctypes.c_char), 
	ctypes.c_double, 
	ctypes.c_double, 
	ctypes.c_double, 
	ctypes.c_double, 
	ctypes.c_double, 
	ctypes.c_double, 
	ctypes.c_double
)
lib.AdjointFTR_getSCVM.restype = np.ctypeslib.ndpointer(dtype=ctypes.c_double, ndim=1, shape=(27,), flags='C')
lib.AdjointFTR_getONmats.argtypes = (
	ctypes.POINTER(ctypes.c_char), 
	ctypes.c_double, 
	ctypes.c_double, 
	ctypes.c_double, 
	ctypes.c_double, 
	ctypes.c_double, 
	ctypes.c_double, 
	ctypes.c_double,
	ctypes.c_double,
	ctypes.c_double
)
lib.AdjointFTR_getONmats.restype = np.ctypeslib.ndpointer(dtype=ctypes.c_double, ndim=1, shape=(12,), flags='C')

class AdjointFTR(ctypes.Structure):
	'''
	Python binding class for c++ source files
	'''
	def __init__(self):
		self.obj = lib.AdjointFTR_new()

	def getSCVM(self, kPerv, Y):
		'''
		Calculate space charge variation matrices given perveance and moments
		'''
		out = lib.AdjointFTR_getSCVM(self.obj, kPerv, Y[0], Y[1], Y[2], Y[3], Y[4], Y[5])

		Mq = out[0:9].reshape((3,3))
		Mp = out[9:18].reshape((3,3))
		Mn = out[18:].reshape((3,3))

		return Mq,Mp,Mn

	def getONmats(self, kPerv, kSol, kQuad, kQuadRot, pipeRadius, Y):
		'''
		Calculate O and N matrices
		'''
		out = lib.AdjointFTR_getONmats(self.obj, kPerv, kSol, kQuad, kQuadRot, pipeRadius, Y[0], Y[1], Y[2], Y[10])

		Omat = out[0:9].reshape((3,3))
		Nmat = out[9:12].reshape((3,1))

		return Omat, Nmat