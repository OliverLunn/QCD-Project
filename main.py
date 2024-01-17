import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

class Kerr_Resonator:

    def __init__(self, n_states, delta, epsilon, K):
        self.n_states = n_states   
        self.delta = delta
        self.epsilon = epsilon
        self.K = K

    def lowering_operator(self):
        n = np.arange(1, self.n_states, 1,dtype=complex)
        operator = np.diag(np.sqrt(n), k=1)
        return operator

    def raising_operator(self):
        n = np.arange(1, self.n_states, 1,dtype=complex)
        operator = np.diag(np.sqrt(n), k=-1)
        return operator

    def n_vector(self, n):
        comp = np.arange(0, self.n_states, 1)
        vector = np.where(comp==n, 1+0j, 0+0j)
        return vector

    def rho_dot_element(self, t, rho, n, m):

        delta = self.delta
        epsilon = self.epsilon
        K = self.K
        i=1j

        rho_d = (-i*delta*n - i*(K/2)*n*(n-1) + i*delta*m + i*(K/2)*m*(m-1) - (1/2)*n - (1/2)*m)*rho[n,m]
        + epsilon*np.sqrt(n)*rho[n-1,m] + epsilon*np.sqrt(m)*rho[n,m-1] - epsilon*np.sqrt(n+1)*rho[n+1,m]
        - epsilon*np.sqrt(m+1)*rho[n,m+1] + np.sqrt(n+1)*np.sqrt(m+1)*rho[n+1,m+1]

        return rho_d
    

    def rho_dot(self, t, rho):
        delta = self.delta
        epsilon = self.epsilon
        K = self.K
        n_states = self.n_states

        rho_dot = np.zeros((n_states, n_states),dtype=complex)
        rho = np.pad(rho, 1, mode="constant",constant_values=0+0j)
        
        for n in range(n_states):
            for m in range(n_states):
                rho_dot[n,m] = self.rho_dot_element(t, rho, n, m)

        return rho_dot
    
    def vectorise(self,mat):
        size = np.shape(mat)[0]
        vec = np.zeros(size**2,dtype=complex)
        for i in range(size):
            for j in range(size):
                vec[j+i*size] = mat[i,j]
        return vec
    
    def unvectorise(self,vec,size):
        mat = np.zeros((size,size),dtype=complex)
        for k in range(size**2):
            i = k//size
            j = k%size
            mat[i,j] = vec[k]
        return mat
    
    def vector_rho_dot(self, t, rho_vec):
        rho = self.unvectorise(rho_vec,self.n_states)
        rho_d = self.rho_dot(t,rho)
        rho_d_vec = self.vectorise(rho_d)
        return rho_d_vec

    def solve(self,rho_0,t_start,t_stop,t_step):
        """
        Perfomring numerical integration of the equations of motion for the density matrix
        """
        soln = integrate.solve_ivp(self.vector_rho_dot,(t_start,t_stop),rho_0,t_eval=np.arange(t_start,t_stop,t_step))
        return soln
    
        
    def photon_number(self, rho, t_idx):
        """
        Calculates the number of photons in the oscillator. 
        """
        a = kr.raising_operator()
        a_dag = kr.lowering_operator()
        expectation_val = np.trace(a_dag@a@rho[:,:,t_idx])

        return expectation_val
    
    def trace_rho_sq(self, rho, t):
        """
        Calculates the trace of the density operator squared
        Takes a 2D array for the denisty operator
        Returns value of trace.
        """
        n_states = self.n_states

        for n in range(n_states):
            for m in range(n_states):
                inner_product = np.inner(rho[:,:]@np.conj(rho[:,:]), self.n_vector(m))
                trace = np.inner(self.n_vector(n), inner_product)

        return trace
    
    
if __name__ == '__main__':

    delta = 0
    epsilon = 1
    K = -0.5
    n_states = 5

    kr = Kerr_Resonator(n_states, delta, epsilon, K)
    rho_0 = np.outer(kr.n_vector(0), kr.n_vector(0))
    rho_0_vec = kr.vectorise(rho_0)

    start, stop = 0, 50
    step=0.01

    soln = kr.solve(rho_0_vec, start, stop, step)
    solution = soln.y
    matrix_solution = np.zeros((n_states,n_states,len(soln.t)), dtype=complex)

    for t in range(len(soln.t)):
        matrix_solution[:,:,t] = kr.unvectorise(solution[:,t], n_states)

    time = np.arange(start, stop, step)
    photon_number = np.zeros(len(time), dtype=complex)   
    tr_rho_sq = np.zeros(len(time), dtype=complex)
    
    for i in range(len(time)):
        tr_rho_sq[i] = kr.trace_rho_sq(matrix_solution[:,:,i], t)
        
    plt.plot(time, tr_rho_sq, ".")
    plt.show()
