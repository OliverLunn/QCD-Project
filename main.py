import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import odeintw

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

        rho_dot = np.zeros((n_states, n_states))
        rho = np.pad(rho, 1, mode="constant")

        for n in range(n_states):
            for m in range(n_states):
                rho_dot[n,m] = self.rho_dot_element(t, rho, n, m)

        return rho_dot


    def solve(self, t_start, t_stop, rho_0):

        soln = odeintw.odeintw(self.rho_dot, rho_0, np.linspace(t_start, t_stop, 100))
        
        return soln
    
if __name__ == '__main__':

    delta = 0
    epsilon = 1
    K = -0.5
    n_states = 3

    kr = Kerr_Resonator(n_states, delta, epsilon, K)
    rho_0 = (np.outer(kr.n_vector(0), kr.n_vector(0))).reshape(-1)

    soln = kr.solve(0, 1, rho_0)
    print(kr.rho_dot(0,rho_0))
    #print(soln)