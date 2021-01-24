# General Markovian two-stage continuous-flow production system w/finite buffer.
# Ref. paper by B. Tan, S.B. Gershwin (2009).
# Building block for a multi-stage production line OR an assembly/disassembly
# line, to be treated with usual decomposition methods.

import numpy as np
import scipy.integrate as spi
import scipy.linalg as spl
from matplotlib import pyplot as plt

class TwoLine():
    """
      Class used to do represent a two-machine line composed by
      two pseudo-machines (Mu, Md) divided by a buffer.
      Also known as "building block" in literature regarding
      decomposition methods.
    """
    
    def set_test_param(self, muu=1.5,mud=1,N=17):
        """
        FOR TESTING. Initializes two-machine line with data from the 
        case study in Tan-Gershwin paper.
        """
        p = 0.005
        r = 0.15
        p1 = 0.015
        r1 = 0.15
        g = 0.01
        h = 0.2
        rQ = 0.15
        
        self.capacity = N

        self.u_TM = np.array([[-g-p, g, p, 0, 0],
                                    [0, -p-h, 0, p, h],
                                    [r, 0, -r, 0, 0],
                                    [0, r, 0, -r, 0],
                                    [rQ, 0, 0, 0, -rQ]])

        self.u_fr = [muu, muu, 0, 0, 0]

        self.d_TM = np.array([[-p1, p1],
                                  [r1, -r1]])

        self.d_fr = [mud, 0]

    def set_param(self, sys, j):
        """
        Initializes two-machine line with data from the original machines
        in the true system.

        Parameters:
        sys    Object representing the original (true) system
        j      Buffer index
        """
        self.capacity = sys.B[j].capacity
        self.u = sys.B[j].upstream
        self.d = sys.B[j].downstream
        self.port = sys.B[j].port
        self.pred, self.succ = sys.B[j].pred, sys.B[j].succ
        self.perc_up = float(sys.B[j].perc_up)
        self.perc_down = float(sys.B[j].perc_down)
        u = self.u-1
        d = self.d-1
        self.u_fr = np.multiply(self.perc_up,sys.M[u].flowrate)
        self.d_fr = np.multiply(self.perc_down,sys.M[d].flowrate)
        self.u_nlocstates = sys.M[u].num_states
        self.d_nlocstates = sys.M[d].num_states
        self.u_states = [*range(sys.M[u].num_states)] # * to unpack range
        self.d_states = [*range(sys.M[d].num_states)]
        self.u_locTM = sys.M[u].TM
        self.d_locTM = sys.M[d].TM

    def class_states(self):
        """
        Classifies pairs of upstream/downstream machine states according to
        their impact on buffer, as per Tan-Gershwin paper.

        incr_state includes states' pairs which increase buffer level;
        decr_state includes states' pairs which decrease buffer level;
        zero_state includes states' pairs which keep buffer level constant.
        """
        tolerance = 1.e-6
        incr_state = []
        m_incr = np.empty(0)
        decr_state = []
        m_decr = np.empty(0)
        zero_state = []
        m_zero = np.empty(0)

        muup = self.u_fr
        mudw = self.d_fr
        
        for j in range(len(mudw)):
            for i in range(len(muup)):
                if muup[i] - mudw[j] > tolerance:
                    incr_state.append([i,j])
                    m_incr = np.append(m_incr,muup[i]-mudw[j])
                elif muup[i] - mudw[j] < -tolerance:
                    decr_state.append([i,j])
                    m_decr = np.append(m_decr,muup[i]-mudw[j])
                else:
                    zero_state.append([i,j])
                    m_zero = np.append(m_zero, 0)

        self.incr_state = incr_state
        self.m_incr = m_incr
        self.decr_state = decr_state
        self.m_decr = m_decr
        self.zero_state = zero_state
        self.m_zero = m_zero

        # set of incr+decr states
        self.incr_decr_state = self.incr_state.copy()
        for i in range(len(self.decr_state)):
            self.incr_decr_state.append(self.decr_state[i])
        M_incr_decr = self.m_incr.copy()
        self.M_incr_decr = np.append(M_incr_decr, self.m_decr)

        # set of decr+zero states, only if both flowrates >0
        self.decr_zero_state = self.decr_state.copy()
        M_decr_zero = self.m_decr.copy()
        for i in range(len(self.zero_state)):
            if self.m_zero[i] > 0:
                self.decr_zero_state.append(self.zero_state[i])
                self.M_decr_zero = np.append(M_decr_zero, self.m_zero[i])

        # set of incr+zero states, only if both flowrates >0
        self.incr_zero_state = self.incr_state.copy()
        M_incr_zero = self.m_incr.copy()
        for i in range(len(self.zero_state)):
            if self.m_zero[i] > 0:
                self.incr_zero_state.append(self.zero_state[i])
                self.M_incr_zero = np.append(M_incr_zero, self.m_zero[i])
                
    def solver(self):
        """
        Main function of TwoLine class.
        Finds the solution of the two-machine line, identifying:
        - Machine state probabilities at steady-state;
        - Average throughput;
        - Average buffer level.

        PLEASE NOTE: errors are expected when machine flowrates tend to be
        close,but not equal.
        ROOT CAUSE: integral of matrix exponential becomes ill-conditioned,
        with complex eigenvalues, therefore generating a garbage solution.
        WORKAROUND: at the moment, only a warning is raised.
        To be verified if the following actions can have positive effect:
        1) decomposing LAMBDA matrix into Jordan form;
        2) increasing numpy float precision.
        """        
        self.create_A1()
        self.create_A2()
        self.create_A3()
        self.create_A4()
        
        self.create_A0()
        self.create_B0()
        self.create_AN()
        self.create_BN()
        print("A1:\n",self.A1)
        print("A2:\n",self.A2)
        print("A3:\n",self.A3)
        print("A4:\n",self.A4)
        print("A0:\n",self.A0)
        print("B0:\n",self.B0)
        print("AN:\n",self.AN)
        print("BN:\n",self.BN)

        OMEGA = -np.dot(spl.inv(self.A4),self.A3)
        LAMBDA = self.A1 - np.linalg.multi_dot([self.A2,spl.inv(self.A4),self.A3])
        #print(np.linalg.cond(LAMBDA))
        #print(np.dot(spl.inv(self.A4),self.A4)) # check identity

        G0 = -np.dot(self.B0, spl.inv(self.A0))
        G0 = G0[:,:len(self.decr_state)]
        E_T0 = -spl.inv(self.A0)
        E_T0 = E_T0[:,:len(self.decr_state)]

        GN = -np.dot(self.BN, spl.inv(self.AN))
        GN = GN[:,:len(self.incr_state)]
        E_TN = -spl.inv(self.AN)
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
    
        eN = spl.expm(LAMBDA*self.capacity)
        #print("eN:\n", eN)
        
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
        inner_state_integral = compute_exp_matrix_integration(LAMBDA, self.capacity)
        #print(np.linalg.cond(inner_state_integral))

        mm = np.concatenate((self.m_incr, self.m_decr))
        eqB = np.dot(mm, inner_state_integral)

        v = (np.ones((1,len(self.incr_decr_state)))+
             np.dot(np.ones((1,len(self.zero_state))),
                    OMEGA))
            
        eqC = c_p0 + c_pN +np.dot(v, inner_state_integral)
        
        lhs = np.vstack((eqA, eqB, eqC))
        rhs = np.zeros((len(self.incr_decr_state),1))
        rhs[-1] = 1

        if 1/np.linalg.cond(lhs) < np.finfo(lhs.dtype).eps:
            w = spl.solve(lhs, rhs)
            print("Condition number is too high!")
##            uu,ss,vv = spl.svd(lhs)
##            ss = np.diag(ss)
##            print(np.shape(uu),np.shape(ss),np.shape(vv))
##            z = np.dot(np.transpose(uu),rhs)
##            y = spl.solve(ss,z)
##            w = np.dot(vv,y)
        else:
            w = spl.solve(lhs, rhs)

        check_empty1 = np.dot(np.concatenate((a,b),axis=1),w)
        check_empty2 = np.linalg.multi_dot([G0,np.concatenate((d,c),axis=1),w])
        check_full1 = np.dot(DECR_a,w)
        check_full2 = np.dot(DECR_b,w)
        print("check_empty working?", np.allclose(check_empty1,check_empty2))
        print("check_full working?", np.allclose(check_full1,check_full2))
        check_w = spl.norm(np.dot(lhs,w)-rhs)

        if check_w > 1.e-5:
            print("Warning! Error of solution is too high:",check_w)
        self.p0 = np.dot(c1_p0,w)
        self.pN = np.dot(c1_pN,w)
        #print("p0:\n", self.p0)
        #print("pN:\n", self.pN)
        #print(self.incr_state,self.pN)
        self.p_inner_incr_decr = np.dot(inner_state_integral,w)
        self.p_inner_zero = np.linalg.multi_dot([OMEGA, inner_state_integral,w])

        self.w = w
        self.LAMBDA = LAMBDA
        self.v = v

        #self.calc_state_prob()
        self.calc_throughput()
        self.calc_exp_buffer()

    def matrix_fill(self, space1, space2, m=[]):
        """
        Calculates sub-matrices A1, A2, A3, A4, according to the input
        parameters.

        Parameters:
        space1      First set of building block states to consider
        space2      Second set of building block states to consider
        """
        size1 = len(space1)
        size2 = len(space2)
        A = np.zeros((size1, size2))
        
        lambda_u = self.u_TM.copy()
        lambda_d = self.d_TM.copy()

        for n in range(size1):
            i, j = space1[n]
            for p in range(size2):
                i1, j1 = space2[p]
                if i1 == i and j1 == j:
                    A[n,p] += lambda_u[i,i]+lambda_d[j,j]
                elif j1 == j:
                    A[n,p] += lambda_u[i1, i]
                elif i1 == i:
                    A[n,p] += lambda_d[j1, j]
                else:
                    A[n,p] += 0
            if len(m) > 0:
                A[n,:] = A[n,:] / m[n]

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

    def matrix_fill_boundary(self, space1, space2, buffer):
        """
        Calculates sub-matrices A0, AN, B0, BN, according to the input
        parameters.

        Parameters:
        space1      First set of building block states to consider
        space2      Second set of building block states to consider
        buffer      "FULL" to generate AN/BN, "EMPTY" to generate A0/B0
        """
        size1 = len(space1)
        size2 = len(space2)
        A = np.zeros((size1, size2))

        lambda_u = self.u_TM.copy()
        lambda_d = self.d_TM.copy()
        
        mu_u = self.u_fr
        mu_d = self.d_fr
        
        if buffer == "EMPTY":
            for n in range(size1):
                i, j = space1[n]
                for p in range(size2):
                    i1, j1 = space2[p]
                    if i1 == i and j1 == j:
                        ld = (mu_u[i]/mu_d[j]) * lambda_d[j,j]
                        A[n,p] += lambda_u[i,i]+ld
                    elif j1 == j:
                        A[n,p] += lambda_u[i1, i]
                    elif i1 == i:
                        ld = (mu_u[i]/mu_d[j1]) * lambda_d[j1,j]
                        A[n,p] += ld
                    else:
                        A[n,p] += 0
            
        if buffer == "FULL":
            for n in range(size1):
                i, j = space1[n]
                for p in range(size2):
                    i1, j1 = space2[p]
                    if i1 == i and j1 == j:
                        lu = (mu_d[j]/mu_u[i]) * lambda_u[i,i]
                        A[n,p] += lu+lambda_d[j,j]
                    elif j1 == j:
                        lu = (mu_d[j]/mu_u[i1]) * lambda_u[i1,i]
                        A[n,p] += lu
                    elif i1 == i:
                        A[n,p] += lambda_d[j1, j]
                    else:
                        A[n,p] += 0

        return A

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
        
    def calc_state_prob(self):
        """
        Calculates machine states' probabilities of two-machine line
        at steady state.
        """
        size_incr = len(self.incr_state)
        size_incr_decr = len(self.incr_decr_state)
        size_zero = len(self.zero_state)
        size_ustates = len(self.u_fr)
        size_dstates = len(self.d_fr)
        self.u_P = np.zeros(size_ustates)
        self.d_P = np.zeros(size_dstates)
        for i in range(size_incr_decr):
            u, d = self.incr_decr_state[i]
            self.u_P[u] += self.p_inner_incr_decr[i]
            self.d_P[d] += self.p_inner_incr_decr[i]
            if i < size_incr:
                self.u_P[u] += self.pN[i]
                self.d_P[d] += self.pN[i]
            else:
                self.u_P[u] += self.p0[i-size_incr]
                self.d_P[d] += self.p0[i-size_incr]
        for i in range(size_zero):
            u, d = self.zero_state[i]
            self.u_P[u] += self.p_inner_zero[i]
            self.d_P[d] += self.p_inner_zero[i]
        
    def calc_prob_empty(self):
        """
        Calculates empty buffer probabilities of two-machine line
        at steady state and records them in array P0.
        In each row, for each combination of upstream state i and
        downstream state j (1st and 2nd column), it reports the
        related empty buffer probability (3rd column).
        """
        size_decr = len(self.decr_state)
        size_ustates = len(self.u_fr)
        size_dstates = len(self.d_fr)
        self.P0 = np.zeros((size_ustates*size_dstates,3))
        cnt=0
        for i in range(size_ustates):
            for j in range(size_dstates):
                self.P0[cnt,0]=i
                self.P0[cnt,1]=j
                for k in range(size_decr):
                    u,d = self.decr_state[k]
                    if u==i and d==j:
                        self.P0[cnt,2]=self.p0[k]
                cnt+=1

    def calc_prob_full(self):
        """
        Calculates full buffer probabilities of two-machine line
        at steady state and records them in array PN.
        In each row, for each combination of upstream state i and
        downstream state j (1st and 2nd column), it reports the
        related full buffer probability (3rd column).
        """
        size_incr = len(self.incr_state)
        size_ustates = len(self.u_fr)
        size_dstates = len(self.d_fr)
        self.PN = self.P0.copy()
        cnt=0
        for i in range(size_ustates):
            for j in range(size_dstates):
                self.PN[cnt,2]=0
                for k in range(size_incr):
                    u,d = self.incr_state[k]
                    if u==i and d==j:
                        self.PN[cnt,2]=self.pN[k]
                cnt+=1
            
    def calc_throughput(self):
        """
        Calculates average throughput of two-machine line at steady state.
        """
        size_incr = len(self.incr_state)
        size_incr_decr = len(self.incr_decr_state)
        size_zero = len(self.zero_state)
        th = 0
        for i in range(size_incr_decr):
            u, d = self.incr_decr_state[i]
            if i < size_incr:
                th += self.d_fr[d]*self.pN[i]
            else:
                th += self.u_fr[u]*self.p0[i-size_incr]
            th += self.u_fr[u]*self.p_inner_incr_decr[i]

        for i in range(size_zero):
            u, d = self.zero_state[i]
            th += self.u_fr[u]*self.p_inner_zero[i]
        
        self.TH = th

    def calc_exp_buffer(self, nbins=100):
        """
        Calculates average buffer of two-machine line
        at steady state.
        """
        N1 = self.capacity*np.sum(self.pN)
        N2 = self.compute_buffer_integration()
        self.exp_buffer = N1+N2

    def compute_buffer_integration(self, nbins=100):
        """
        Integral formula for calculating average buffer level.
        Comparable to "compute_exp_matrix_integration" function.
        """
        f = lambda x: x*np.linalg.multi_dot([self.v,
                                             spl.expm(self.LAMBDA*x),
                                             self.w])
        xv = np.linspace(0,self.capacity,nbins)
        result = np.apply_along_axis(f,0,xv.reshape(1,-1))
        N = np.trapz(result,xv)
        return N
    
def compute_exp_matrix_integration(A,T,nbins=100):
    """
    Integral formula for calculating integral of matrix exponential.

    Kudos to user ArtificiallyIntelligence on Mathematics Stack Exchange.
    https://math.stackexchange.com/questions/658276/integral-of-matrix-exponential
    """
    f = lambda x: spl.expm(A*x)
    xv = np.linspace(0,T,nbins)
    result = np.apply_along_axis(f,0,xv.reshape(1,-1))
    return np.trapz(result,xv)

## initialize inputs
def main():
    up = 0.9
    low = 0.9
    linea = TwoLine()
    muu_space = np.linspace(low,up,round((up-low)/0.1)+1)
    E = np.empty_like(muu_space)
    th = np.empty_like(muu_space)
    for i in range(0,len(muu_space)):
        m = muu_space[i]
        linea.set_test_param(muu=m)
        linea.class_states()
        linea.solver()
        E[i] = linea.exp_buffer
        th[i] = linea.TH
        print(m, E[i], th[i])

    if len(muu_space) > 1:
        plt.plot(muu_space,E)
        plt.plot(muu_space,th)
        plt.show()
        #print("Tempo calcolo (s):", e-s)

if __name__ == "__main__":
    main()
