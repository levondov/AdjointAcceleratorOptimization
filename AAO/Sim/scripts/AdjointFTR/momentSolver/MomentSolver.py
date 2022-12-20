import sys, os, pathlib
from pathlib import Path
import numpy as np
simRoot = Path(__file__).parent.parent
sys.path.append( os.path.join(simRoot,'pythonBindings') )
print(simRoot)

from adjointFTRBindings import *
import numpy as np
import copy 

class MomentSolver:

    def __init__(self, lattice=None, energy=5e3, current=0.0, pipeRadius=0.0, zInterval=(0.0, 0.322), stepSize=0.0001):
        self.initialMoments = self.GetInitialMoments()
        self.energy = energy
        self.current = current
        self.pipeRadius = pipeRadius
        self.lattice = copy.deepcopy(lattice)
        self.latticeDefault = copy.deepcopy(lattice)

        self.zInterval = zInterval
        self.stepSize = stepSize

        self.MomentSolverUtility = MomentSolverUtility()
        self.rho, self.kPerv = self.MomentSolverUtility.GetParameters(self.energy, self.current)

        # c++ bindings
        self.bindings = AdjointFTR()

    def GetInitialMoments(self):
        '''
        return starting conditions for the moment equations (based off Santiago's values)
        can be updated whenever.
        '''
        Q_plus = 0.5*(2.2581**2*1e-6 + 0.2258**2*1e-6)
        Q_minus = 0.5*(2.2581**2*1e-6 - 0.2258**2*1e-6)
        Q_x = 0.0
        P_plus = 0
        P_minus = 0
        P_x = 0
        E_plus = (7.0855**2*1e-6 + 0.70855**2*1e-6)
        E_minus = (7.0855**2*1e-6 - 0.70855**2*1e-6)
        E_x = 0
        L = 0
        phi = 0
        return np.array([Q_plus,Q_minus,Q_x,P_plus,P_minus,P_x,E_plus,E_minus,E_x,L,phi])

    def UpdateLattice(self, lattice=None, params=None):
        '''
        Update lattice parameters
        '''

        if ( type(lattice) is np.ndarray ):
            self.lattice = copy.deepcopy(lattice)
            self.latticeDefault = copy.deepcopy(lattice)

        self.optParams = params

        if ( type(params) is np.ndarray ):
            for i,elem in enumerate(self.lattice):
                self.lattice[i]['zstart'] = self.latticeDefault[i]['zstart'] * params[i]['zstart']
                self.lattice[i]['zend'] = self.lattice[i]['zstart'] + self.latticeDefault[i]['length']
                self.lattice[i]['rotation'] = self.latticeDefault[i]['rotation'] * params[i]['rotation']
                self.lattice[i]['dbdx'].SetScaleFactor( params[i]['dbdx'] )

    def Run(self, verbose=False):
        '''
        Run the moment equations over a given interval with a given step size
        '''

        # run forward integration of moment equations
        z = np.arange(self.zInterval[0], self.zInterval[1] + self.stepSize, self.stepSize) # all steps

        odefunc = lambda z,Y : self.OdeMoments(z,Y,self.lattice)
        y,ksol,kquad,Omat,Nmat = self.MomentSolverUtility.ode3(odefunc,self.zInterval[0], self.stepSize, self.zInterval[1], self.initialMoments, verbose=verbose)

        # grab solutions
        self.z = z
        self.y = y
        self.ksol = ksol
        self.kquad = kquad
        self.Omat = Omat
        self.Nmat = Nmat
        
        return z,y,ksol,kquad

    def RunAdjoint(self, verbose=False):
        '''
        Run the adjoint moment equations backwards
        '''

        # run forward integration of moment equations
        z = np.arange(self.zInterval[1],self.zInterval[0] - self.stepSize,-self.stepSize) # all steps

        # setup initial conditions
        _, _, self.initialMomentsAdjoint = self.GetFoM_And_DFoM()
        initialMoments = np.concatenate((self.initialMomentsAdjoint, np.array([self.y[-1,-1]]), self.y[:,-1]))

        odefunc = lambda z,Y : self.OdeMomentsAdjoint(z,Y,self.lattice)
        y,ksol,kquad,_,_ = self.MomentSolverUtility.ode3(odefunc,self.zInterval[1], -self.stepSize, self.zInterval[0], initialMoments, verbose=verbose)

        self.zAdj = np.flip(z)
        self.yAdj = np.flip(y, 1)

        return self.zAdj, self.yAdj, ksol, kquad

    def OdeMoments(self, z, Y, lattice):
        '''
        Main function that solves Tom's moment equations
        
        Y(1) - Q+
        Y(2) - Q-
        Y(3) - Qx
        Y(4) - P+
        Y(5) - P-
        Y(6) - Px
        Y(7) - E+
        Y(8) - E-
        Y(9) - Ex
        Y(10) - L
        '''

        k_perv = self.kPerv
        r_pipe = self.pipeRadius
        rho = self.rho

        k_sol, k_quad, psi = self.GetMagnetFields(z, lattice, rho)
                
        # use bindings to calculate O and N matrices
        O_mat,N_mat = self.bindings.getONmats(k_perv, k_sol, k_quad, psi, r_pipe, Y)
        
        # System of 10 equations
        dydt = np.array([
        # dQ/dz
        Y[3],
        Y[4],
        Y[5],
        # dP/dz
        Y[6] + np.matmul( O_mat[0,:], np.reshape(Y[0:3],(3,1)) )[0],
        Y[7] + np.matmul( O_mat[1,:], np.reshape(Y[0:3],(3,1)) )[0],
        Y[8] + np.matmul( O_mat[2,:], np.reshape(Y[0:3],(3,1)) )[0],
        # dE/dz
        np.matmul( O_mat[0,:], np.reshape(Y[3:6],(3,1)) )[0] + N_mat[0,0]*Y[9],
        np.matmul( O_mat[1,:], np.reshape(Y[3:6],(3,1)) )[0] + N_mat[1,0]*Y[9],
        np.matmul( O_mat[2,:], np.reshape(Y[3:6],(3,1)) )[0] + N_mat[2,0]*Y[9],
        # dL/dz
        -1.0*np.matmul( N_mat.T, Y[0:3] )[0],
        # dphi/dz
        -1.0*k_sol/2.0
        ])
        
        return dydt,k_sol,k_quad,O_mat,N_mat

    def OdeMomentsAdjoint(self, z, Yt, lattice):  
        '''
        Main function to solve the adjoint equations (backwards)
        
        Here we have the original 11 moments + the 11 adjoint moments for a total of a 22 variable ode solve
        '''
        Y = Yt[0:11] # adjoint variables
        Y2 = Yt[11:] # moment variables
        
        k_perv = self.kPerv
        r_pipe = self.pipeRadius
        rho = self.rho
        
        k_sol, k_quad, psi = self.GetMagnetFields(z, lattice, rho)
                
        # use bindings to calculate O and N matrices
        O_mat,N_mat = self.bindings.getONmats(k_perv, k_sol, k_quad, psi, r_pipe, Y2)
        
        # Calculate special matrix due to space charge variations
        Mq,Mp,Mn = self.bindings.getSCVM(k_perv, Y2)

        # adjoint dE_dot / dz special
        dedot = np.matmul( Y[3:6], Mq ) + Y[9] * np.matmul( Y2[0:3], Mn ) - np.matmul( Y[0:3], Mp ) - np.matmul( Y[0:3], Mn ) * Y2[9]

        # System of 20 equations
        dydt = np.array([
        # adjoint dQ/dz
        Y[3],
        Y[4],
        Y[5],
        # adjoint dP/dz
        Y[6] + np.matmul( O_mat[0,:], np.reshape(Y[0:3],(3,1)) )[0],
        Y[7] + np.matmul( O_mat[1,:], np.reshape(Y[0:3],(3,1)) )[0],
        Y[8] + np.matmul( O_mat[2,:], np.reshape(Y[0:3],(3,1)) )[0],
        # adjoint dE/dz
        np.matmul( O_mat[0,:], np.reshape(Y[3:6],(3,1)) )[0] + N_mat[0,0]*Y[9] + dedot[0],
        np.matmul( O_mat[1,:], np.reshape(Y[3:6],(3,1)) )[0] + N_mat[1,0]*Y[9] + dedot[1],
        np.matmul( O_mat[2,:], np.reshape(Y[3:6],(3,1)) )[0] + N_mat[2,0]*Y[9] + dedot[2],
        # adjoint dL/dz
        -1.0*np.matmul( N_mat.T, Y[0:3] )[0],        
        # adjoint dphi/dz
        -1.0*k_sol/2.0,        
        # dQ/dz
        Y2[3],
        Y2[4],
        Y2[5],
        # dP/dz
        Y2[6] + np.matmul( O_mat[0,:], np.reshape(Y2[0:3],(3,1)) )[0],
        Y2[7] + np.matmul( O_mat[1,:], np.reshape(Y2[0:3],(3,1)) )[0],
        Y2[8] + np.matmul( O_mat[2,:], np.reshape(Y2[0:3],(3,1)) )[0],
        # dE/dz
        np.matmul( O_mat[0,:], np.reshape(Y2[3:6],(3,1)) )[0] + N_mat[0,0]*Y2[9],
        np.matmul( O_mat[1,:], np.reshape(Y2[3:6],(3,1)) )[0] + N_mat[1,0]*Y2[9],
        np.matmul( O_mat[2,:], np.reshape(Y2[3:6],(3,1)) )[0] + N_mat[2,0]*Y2[9],
        # dL/dz
        -1.0*np.matmul( N_mat.T, Y2[0:3] )[0],
        # dphi/dz
        -1.0*k_sol/2.0,
        ])

        return dydt,k_sol,k_quad,O_mat,N_mat        

    def GetCOM(self, y):
        '''
        Calculate the constant of motion
        '''
        L = y[9,:]
        EQ = y[6,:]*y[0,:] + y[7,:]*y[1,:] + y[8,:]*y[2,:]
        PP = y[3,:]**2 + y[4,:]**2 + y[5,:]**2
        motion = EQ + (0.5)*L**2 - (0.5)*PP
        
        return motion

    def GetFoM_And_DFoM(self, index = -1, k0 = 7):
        '''
        calculates the Figure of Merit + adjoint equation initial conditions given a set of moment values
        
        i.e.
        Given Q+,Q-,Qx,P+,P-,Px,E+,E-,Ex,L , calculate FOM:
        '''        
        komega = self.ksol[index]
        y = self.y
        e1 = 1.0
        e2 = 0.0

        f5_tmp = y[6,index] + 0.5 * komega**2 * y[0,index] - komega * y[9,index]
        f4_tmp = y[6,index] - 0.5 * komega**2 * y[0,index] + self.kPerv      

        # figure of merit broken into pieces for ease of reading
        FoM1 = 0.5 * ( y[4,index]**2 + y[5,index]**2 + y[3,index]**2 )
        FoM2 = 0.5 * (k0**2) * ( y[1,index]**2 + y[2,index]**2 )
        FoM3 = 0.5 * (k0**(-2)) * ( y[7,index]**2 + y[8,index]**2 )
        FoM4 = 0.5 * (k0**(-2)) * e1 * f4_tmp**2
        FoM5 = 0.5 * (k0**(-2)) * e2 * f5_tmp**2
            
        FoM = FoM1 + FoM2 + FoM3 + FoM4 + FoM5
        FoMp = np.array([FoM1,FoM2,FoM3,FoM4,FoM5])  

        # now calculate derivatives for adjoint variables
        # adjoint variables calculated from FoM
        dP_p = y[3,index]
        dP_m = y[4,index]
        dP_x = y[5,index]
        
        dE_p = (k0**(-2)) * e2 * (f5_tmp) * (0.5*komega**2) * (-1) + (k0**(-2)) * e1 * (f4_tmp) * (-0.5*komega**2) * (-1)
        dE_m = -k0**(2)*y[1,index]
        dE_x = -k0**(2)*y[2,index]
        
        dQ_p = (k0**(-2)) * e2 * (f5_tmp) * (-1) + (k0**(-2)) * e1 * (f4_tmp) * (-1)
        dQ_m = -k0**(-2)*y[7,index]
        dQ_x = -k0**(-2)*y[8,index]
        
        dL = (k0**(-2)) * e2 * (f5_tmp) * (-1) * komega
        
        dFoM = np.array([dQ_p,dQ_m,dQ_x,dP_p,dP_m,dP_x,dE_p,dE_m,dE_x,dL])          
            
        return FoM, FoMp, dFoM

    def GetDF(self):
        '''
        Calculate Omat,Nmat perturbation matrices
        '''
        pert = 0.01
        OmatDefault = self.Omat
        NmatDefault = self.Nmat
        paramsDefault = copy.deepcopy(self.optParams)
        params = None
        df = []
        for i,elem in enumerate(paramsDefault):
            for elemName in elem:
                params = copy.deepcopy(paramsDefault)
                params[i][elemName] = paramsDefault[i][elemName] + paramsDefault[i][elemName] * pert
                self.UpdateLattice(params = params)
                Omat,Nmat = self.CalcON(self.lattice)
                for ii in range(len(self.z)):
                    Omat[ii] = Omat[ii] - OmatDefault[ii]
                    Nmat[ii] = Nmat[ii] - NmatDefault[ii]  

                tmp = self.CalcInt(Omat, Nmat)
                df.append(tmp)

        return np.array(df)

    def CalcON(self, lattice):
        '''
        Calculate O and N matrices
        '''
        O_mat = np.empty(len(self.z), dtype=np.ndarray)
        N_mat = np.empty(len(self.z), dtype=np.ndarray)

        k_perv = self.kPerv
        r_pipe = self.pipeRadius
        rho = self.rho

        for i,z in enumerate(self.z):
            Y = self.y[:,i]
            
            k_sol, k_quad, psi = self.GetMagnetFields(z, lattice, rho)
        
            # use bindings to calculate O and N matrices
            O_mat[i],N_mat[i] = self.bindings.getONmats(k_perv, k_sol, k_quad, psi, r_pipe, Y)

        return O_mat, N_mat

    def CalcInt(self, Omat, Nmat):
        '''
        Calculate integral
        '''
        int1 = np.zeros(len(self.z))
        int2 = np.zeros(len(self.z))
        int3 = np.zeros(len(self.z))
        int4 = np.zeros(len(self.z))

        y = self.y
        yadj = self.yAdj

        for i in range(len(self.z)):
            int1[i] = np.matmul( yadj[3:6,i], np.matmul(Omat[i], y[0:3,i]) )
            int2[i] = yadj[9,i] * np.matmul( y[0:3,i], Nmat[i] )
            int3[i] = -1 * np.matmul( yadj[0:3,i], np.matmul(Omat[i], y[3:6,i]) )
            int4[i] = -1 * np.matmul( yadj[0:3,i], Nmat[i] * y[9,i] )

        return np.trapz(int1+int2+int3+int4, self.z)

    def GetMagnetFields(self, z, lattice, rho):

        # are we in magnets? get strengths
        k_sol = 0.0
        k_quad = 0.0
        psi = 0.0        
        for elem in lattice:
            # find solenoid
            if ( elem['type'] == 'solenoid' ):
                # are we within the field?
                if ( z >= elem['zstart'] and z < elem['zend'] ):
                    # add the field
                    k_sol = k_sol + (elem['dbdx'].GetValue( z - elem['zstart'] ) / rho)

            # find quadrupole
            if ( elem['type'] == 'quad' ):
                # are we within the field?
                if ( z >= elem['zstart'] and z < elem['zend'] ):
                    # add the field
                    k_quad = k_quad + (elem['dbdx'].GetValue( z - elem['zstart'] ) / rho)
                    psi = elem['rotation']

        return k_sol, k_quad, psi

class MomentSolverUtility:

    def __init__(self):
        pass

    def GetParameters(self, energy=5e3,current=0.0):
        '''
        Grab beam ridgidity
        '''
        # Parameters
        e         = 1.60217733E-19 #C
        m         = 9.1093897E-31 #kg
        Energy    = energy # eV
        c         = 2.997924E8 # m/s

        gamma     = 1+((Energy)/(510998.9461));
        beta      = np.sqrt((gamma*gamma)-1)/gamma
        v         = beta*c
        bg        = beta*gamma
        rho       = bg*c*(m/e) 
        
        k_perv = (1.0/(4.0*np.pi))*(c*377.0*current) / (m*v**3*gamma**3/e);   
        
        return rho,k_perv

    def ode3(self,F,t0,h,tfinal,y0,verbose=False):
        '''
        third order classical Runge-Kutta ODE solver
        '''
        y = y0
        tsteps = np.arange(t0, tfinal, h)           

        # for extra params
        ksol = np.zeros(len(tsteps)+1)
        kquad = np.zeros(len(tsteps)+1)
        Omat = np.empty(len(tsteps)+1, dtype=np.ndarray)
        Nmat = np.empty(len(tsteps)+1, dtype=np.ndarray)

        _,sol,quad,Om,Nm = F(t0,y0)
        ksol[0] = sol
        kquad[0] = quad
        Omat[0] = Om
        Nmat[0] = Nm

        yout = np.zeros((len(y0),len(tsteps)+1))
        yout[:,0] = y

        if verbose:
            N = len(tsteps)
            NN = int(N/10.0)
            for i,t in enumerate(tsteps):
                t1,sol,quad,Om,Nm = F(t,y)
                s1 = h*t1
                t2,_,_,_,_ = F(t+h/2.0, y+s1/2.0)
                s2 = h*t2
                t3,_,_,_,_ = F(t+h, y-s1+2.0*s2)
                s3 = h*t3
                y = y + (s1 + 4.0*s2 + s3)/6.0
                yout[:,i+1] = y
                ksol[i+1] = sol
                kquad[i+1] = quad
                Omat[i+1] = Om
                Nmat[i+1] = Nm
                if i % NN == 0:
                    zstr = ' z = ' + str(round(t,4))[0:6]
                    print('{:<15}'.format(zstr) +  ' | ' + str(int(1.0*i/NN*10)) + '%')

            zstr = ' z = ' + str(round(t+h,4))[0:6]
            print('{:<15}'.format(zstr) +  ' | 100%')            
        else:
            for i,t in enumerate(tsteps):
                t1,sol,quad,Om,Nm = F(t,y)
                s1 = h*t1
                t2,_,_,_,_ = F(t+h/2.0, y+s1/2.0)
                s2 = h*t2
                t3,_,_,_,_ = F(t+h, y-s1+2.0*s2)
                s3 = h*t3
                y = y + (s1 + 4.0*s2 + s3)/6.0
                yout[:,i+1] = y
                ksol[i+1] = sol
                kquad[i+1] = quad
                Omat[i+1] = Om
                Nmat[i+1] = Nm             
                
        return yout,ksol,kquad,Omat,Nmat