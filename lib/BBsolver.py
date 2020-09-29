# General Markovian two-stage continuous-flow production system w/finite buffer.
# Ref. paper by B. Tan, S.B. Gershwin (2009).
# Building block for a multi-stage production line OR an assembly/disassembly
# line, to be treated with usual decomposition methods.
# 
# N.B.: errors are expected when mu_u and mu_d tend to be close (but not equal).
# CAUSE: integral of matrix exponential becoming ill-conditioned, with complex
# eigenvalues, therefore generating garbage solution.
# SOLUTION: at the moment, only a warning is raised.
# To be verified if decomposing LAMBDA into Jordan form can have positive effect
# https://math.stackexchange.com/questions/658276/integral-of-matrix-exponential

import numpy as np
import scipy.integrate as spi
import scipy.linalg as spl
import pandas as pd
import time
from matplotlib import pyplot as plt

class BB:
    def set_param(self, muu=1.5,mud=5.5,N=17):
        p = 0.005
        r = 0.15
        p1 = 0.015
        r1 = 0.15
        g = 0.01
        h = 0.2
        rQ = 0.15
        
        self.buffer = N

        self.statename_u = ["1","-1","D1","D-1","DQ"]
        self.statename_d = ["1'","0'"]

        self.lambda_u = np.array([[-g-p, g, p, 0, 0],
                                    [0, -p-h, 0, p, h],
                                    [r, 0, -r, 0, 0],
                                    [0, r, 0, -r, 0],
                                    [rQ, 0, 0, 0, -rQ]])

        self.mu_u = [muu, muu, 0, 0, 0]

        self.mu_d = [mud, 0]

        self.lambda_d = np.array([[-p1, p1],
                                  [r1, -r1]])
    

    def class_states(self):
        tolerance = 0.000001
        incr_state = []
        m_incr = np.empty(0)
        
        decr_state = []
        m_decr = np.empty(0)
        
        zero_state = []
        m_zero = np.empty(0)

        muup = self.mu_u
        mudw = self.mu_d
        stateup = self.statename_u
        statedw = self.statename_d
        
        for j in range(len(mudw)):
            for i in range(len(muup)):
                if muup[i] - mudw[j] > tolerance:
                    #incr_state.append([stateup[i],statedw[j]])
                    incr_state.append([i,j])
                    m_incr = np.append(m_incr,muup[i]-mudw[j])
                elif muup[i] - mudw[j] < -tolerance:
                    #decr_state.append([stateup[i],statedw[j]])
                    decr_state.append([i,j])
                    m_decr = np.append(m_decr,muup[i]-mudw[j])
                else:
                    #zero_state.append([stateup[i],statedw[j]])
                    zero_state.append([i,j])
                    m_zero = np.append(m_zero, 0)

        self.incr_state = incr_state
        self.m_incr = m_incr
        self.decr_state = decr_state
        self.m_decr = m_decr
        self.zero_state = zero_state
        self.m_zero = m_zero

        # insieme stati incr+decr
        self.incr_decr_state = self.incr_state.copy()
        for i in range(len(self.decr_state)):
            self.incr_decr_state.append(self.decr_state[i])
        M_incr_decr = self.m_incr.copy()
        self.M_incr_decr = np.append(M_incr_decr, self.m_decr)

        # insieme stati decr+zero (S0), solo se entrambi mu != 0
        self.decr_zero_state = self.decr_state.copy()
        M_decr_zero = self.m_decr.copy()
        for i in range(len(self.zero_state)):
            if self.m_zero[i] > 0:
                self.decr_zero_state.append(self.zero_state[i])
                self.M_decr_zero = np.append(M_decr_zero, self.m_zero[i])

        self.incr_zero_state = self.incr_state.copy()
        M_incr_zero = self.m_incr.copy()
        for i in range(len(self.zero_state)):
            if self.m_zero[i] > 0:
                self.incr_zero_state.append(self.zero_state[i])
                self.M_incr_zero = np.append(M_incr_zero, self.m_zero[i])


    def matrix_fill(self, space1, space2, m=[]):
        size1 = len(space1)
        size2 = len(space2)
        A = np.zeros((size1, size2))
        
        lambda_u = self.lambda_u
        lambda_d = self.lambda_d

        for n in range(size1):
            i, j = space1[n]
            for p in range(size2):
                i1, j1 = space2[p]
                if i1 == i and j1 == j:
                    A[n,p] = lambda_u[i,i]+lambda_d[j,j]
                elif j1 == j:
                    A[n,p] = lambda_u[i1, i]
                elif i1 == i:
                    A[n,p] = lambda_d[j1, j]
                else:
                    A[n,p] = 0
            if len(m) > 0:
                A[n,:] = A[n,:] / m[n]

        return A

    def matrix_fill_boundary(self, space1, space2, buffer):
        size1 = len(space1)
        size2 = len(space2)
        A = np.zeros((size1, size2))

        lambda_u = self.lambda_u.copy()
        lambda_d = self.lambda_d.copy()
        
        mu_u = self.mu_u
        mu_d = self.mu_d
        
        if buffer == "EMPTY":
            for n in range(size1):
                i, j = space1[n]
                for p in range(size2):
                    i1, j1 = space2[p]
                    if i1 == i and j1 == j:
                        lambda_d[j,j] = (mu_u[i]/mu_d[j]) * lambda_d[j,j]
                        A[n,p] = lambda_u[i,i]+lambda_d[j,j]
                    elif j1 == j:
                        A[n,p] = lambda_u[i1, i]
                    elif i1 == i:
                        lambda_d[j1,j] = (mu_u[i]/mu_d[j1]) * lambda_d[j1,j]
                        A[n,p] = lambda_d[j1, j]
                    else:
                        A[n,p] = 0
            
        if buffer == "FULL":
            for n in range(size1):
                i, j = space1[n]
                for p in range(size2):
                    i1, j1 = space2[p]
                    if i1 == i and j1 == j:
                        lambda_u[i,i] = (mu_d[j]/mu_u[i]) * lambda_u[i,i]
                        A[n,p] = lambda_u[i,i]+lambda_d[j,j]
                    elif j1 == j:
                        lambda_u[i1,i] = (mu_d[j]/mu_u[i1]) * lambda_u[i1,i]
                        A[n,p] = lambda_u[i1, i]
                    elif i1 == i:
                        A[n,p] = lambda_d[j1, j]
                    else:
                        A[n,p] = 0

        return A

    def create_A1(self):
        incr_decr_state = self.incr_decr_state
        
        self.A1 = self.matrix_fill(incr_decr_state,incr_decr_state, self.M_incr_decr)
        
    def create_A2(self):
        incr_decr_state = self.incr_decr_state
        zero_state = self.zero_state

        self.A2 = self.matrix_fill(incr_decr_state,zero_state,self.M_incr_decr)

    def create_A3(self):
        incr_decr_state = self.incr_decr_state
        zero_state = self.zero_state

        self.A3 = self.matrix_fill(zero_state, incr_decr_state)

    def create_A4(self):
        zero_state = self.zero_state

        self.A4 = self.matrix_fill(zero_state, zero_state)

    def create_A0(self):
        decr_zero_state = self.decr_zero_state

        self.A0 = self.matrix_fill_boundary(decr_zero_state,
                                            decr_zero_state, "EMPTY")

    def create_B0(self):
        incr_state = self.incr_state
        decr_zero_state = self.decr_zero_state

        self.B0 = self.matrix_fill_boundary(incr_state,
                                            decr_zero_state, "EMPTY")

    def create_AN(self):
        incr_zero_state = self.incr_zero_state

        self.AN = self.matrix_fill_boundary(incr_zero_state,
                                            incr_zero_state, "FULL")

    def create_BN(self):
        decr_state = self.decr_state
        incr_zero_state = self.incr_zero_state

        self.BN = self.matrix_fill_boundary(decr_state,
                                            incr_zero_state, "FULL")

    def calc_throughput(self):
        size_incr = len(self.incr_state)
        size_incr_decr = len(self.incr_decr_state)
        size_zero = len(self.zero_state)
        th = 0
        for i in range(size_incr_decr):
            u, d = self.incr_decr_state[i]
            if i < size_incr:
                th = th + self.mu_d[d]*self.pN[i]
            else:
                th = th +self.mu_u[u]*self.p0[i-size_incr]
            th = th + self.mu_u[u]*self.p_inner_incr_decr[i]

        for i in range(size_zero):
            u, d = self.zero_state[i]
            th = th + self.mu_u[u]*self.p_inner_zero[i]
        
        self.TH = th

    def compute_buffer_integration(self, nbins=100):
        f = lambda x: x*np.linalg.multi_dot([self.v,
                                             spl.expm(self.LAMBDA*x),
                                             self.w])
        xv = np.linspace(0,self.buffer,nbins)
        result = np.apply_along_axis(f,0,xv.reshape(1,-1))
        N = np.trapz(result,xv)
        return N

    def calc_exp_buffer(self, nbins=100):
        N1 = self.buffer*np.sum(self.pN)
        N2 = self.compute_buffer_integration()
        self.exp_buffer = N1+N2
    
    def solver(self):
        self.class_states()
        
        self.create_A1()
        self.create_A2()
        self.create_A3()
        self.create_A4()
        
        self.create_A0()
        self.create_B0()
        self.create_AN()
        self.create_BN()

        LAMBDA = self.A1 - np.linalg.multi_dot([self.A2,np.linalg.inv(self.A4),self.A3])
        OMEGA = -np.dot(np.linalg.inv(self.A4),self.A3)
        D, U = np.linalg.eig(LAMBDA)

        G0 = -np.dot(self.B0, np.linalg.inv(self.A0))
        G0 = G0[:,:len(self.decr_state)]
        E_T0 = -np.linalg.inv(self.A0)
        E_T0 = E_T0[:,:len(self.decr_state)]

        GN = -np.dot(self.BN, np.linalg.inv(self.AN))
        GN = GN[:,:len(self.incr_state)]
        E_TN = -np.linalg.inv(self.AN)
        E_TN = E_TN[:,:len(self.incr_state)]

        a = np.diag(self.m_incr)
        b = np.zeros((len(self.incr_state),len(self.decr_state)))
        c = -np.diag(self.m_decr)
        d = np.zeros((len(self.decr_state),len(self.incr_state)))
        INCR = np.concatenate((a,b),axis=1) - np.dot(G0, np.concatenate((d,c),axis=1))

        
        INCR = INCR[1:INCR.shape[0],:]
        c1_p0 = np.linalg.multi_dot([E_T0,
                                    np.concatenate((d,c),axis=1)])
        c_p0 = np.linalg.multi_dot([np.ones((1,E_T0.shape[0])),
                                    E_T0,
                                    np.concatenate((d,c),axis=1)])
    
        eN = spl.expm(LAMBDA*self.buffer)
        
        DECR_a =  np.dot(np.concatenate((d,c),axis=1),eN)
        DECR_b = np.dot(GN, np.concatenate((a,b),axis=1))
        DECR_b = np.dot(DECR_b, eN)
        DECR = DECR_a - DECR_b
        
        DECR = DECR[:DECR.shape[0]-1,:]
        c1_pN = np.linalg.multi_dot([E_TN,
                                    np.concatenate((a,b),axis=1),
                                    eN])
        c_pN = np.linalg.multi_dot([np.ones((1,E_TN.shape[0])),
                                    E_TN,
                                    np.concatenate((a,b),axis=1),
                                    eN])
        
        eqA = np.vstack((INCR,DECR))
    
        #THIS DID NOT WORK (generated negative probabilities due to ill-conditioned eN)
        #inner_state_integral = np.dot(np.linalg.inv(LAMBDA),
        #                              (eN-np.identity(len(LAMBDA))))
        
        #THIS FUNCTION SEEMS TO BE THE PROPER VERSION!
        #Kudos to ArtificiallyIntelligence on Mathematics Stack Exchange.
        #https://math.stackexchange.com/questions/658276/integral-of-matrix-exponential
        inner_state_integral = compute_exp_matrix_integration(LAMBDA, self.buffer)
        print(np.linalg.cond(inner_state_integral))

        mm = np.concatenate((self.m_incr, self.m_decr))
        eqB = np.dot(mm, inner_state_integral)    

        v = (np.ones((1,len(self.incr_decr_state)))+
             np.dot(np.ones((1,len(self.zero_state))),
                    OMEGA))
            
        eqC = c_p0 + c_pN +np.dot(v, inner_state_integral)
        
        left_side = np.vstack((eqA, eqB, eqC))
        right_side = np.zeros((len(self.incr_decr_state),1))
        right_side[len(self.incr_decr_state)-1] = 1        

        #w = np.dot(np.linalg.inv(left_side), right_side)
        w = np.linalg.solve(left_side, right_side)

        check_w = np.linalg.norm(np.dot(left_side,w)-right_side)
        
        if check_w > 1.e-6:
            print("Warning! Error of solution is too high:",
              check_w)
        print(check_w)

        self.p0 = np.dot(c1_p0,w)
        self.pN = np.dot(c1_pN,w)
        self.p_inner_incr_decr = np.dot(inner_state_integral,w)
        self.p_inner_zero = np.linalg.multi_dot([OMEGA, inner_state_integral,w])

        self.w = w
        self.LAMBDA = LAMBDA
        self.v = v

        self.calc_throughput()
        self.calc_exp_buffer()
    
def compute_exp_matrix_integration(A,T,nbins=100):
    f = lambda x: spl.expm(A*x)
    xv = np.linspace(0,T,nbins)
    result = np.apply_along_axis(f,0,xv.reshape(1,-1))
    return np.trapz(result,xv)

## initialize inputs
low= 5
up = 5
linea = BB()
muu_space = np.linspace(low,up,round((up-low)/0.1)+1)
E = np.empty_like(muu_space)
th = np.empty_like(muu_space)
for i in range(0,len(muu_space)):
    m = muu_space[i]
    s = time.time()
    linea.set_param(muu=m)
    linea.solver()
    E[i] = linea.exp_buffer
    th[i] = linea.TH
    print(m, E[i], th[i])
    e = time.time()

if len(muu_space) > 1:
    plt.plot(muu_space,E)
    plt.plot(muu_space,th)
    plt.show()
    #print("Tempo calcolo (s):", e-s)
