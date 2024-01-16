import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import odeintw

class Kerr_Resonator:

    def __init__(self, n_states,delta,epsilon,K):

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
        comp = np.arange(0,self.n_states,1)
        vector = np.where(comp==n,1+0j,0+0j)
        return vector

    def rho_dot(self, t, rho):
        delta = self.delta
        epsilon = self.epsilon
        K = self.K

        a = self.lowering_operator()
        a_dag = self.raising_operator()
        i = 1j #imaginary unit

        rho_d = - i*delta*np.dot(a@a_dag@a,rho) - i*(K/2)*np.dot(a@a_dag@a_dag@a@a,rho) - epsilon*np.dot(a@a_dag,rho)
        + epsilon*np.dot(a@a,rho)+ i*delta*np.dot(np.dot(a,rho),a_dag@a) + i*(K/2)*np.dot(np.dot(a,rho),a_dag@a_dag@a@a) - epsilon*np.dot(np.dot(a,rho),a_dag)
        + epsilon*np.dot(np.dot(a,rho),a) - (1/2)*(np.dot(a_dag@a,rho) - 2*np.dot(a,rho)@a_dag) + np.dot(rho,a_dag@a)
    
        return rho_d

    
    def rho_dot_element(self,t,rho,n,m):
        delta = self.delta
        epsilon = self.epsilon
        K = self.K
        i=1j
        return (-i*delta*n - i*(K/2)*n(n-1) + i*delta*m + i*(K/2)*m(m-1) - (1/2)*n - (1/2)*m)*rho[n,m] + epsilon*np.sqrt(n)*rho[n-1,m] + epsilon*np.sqrt(m)*rho[n,m-1] - epsilon*np.sqrt(n+1)*rho[n+1,m] - epsilon*np.sqrt(m+1)*rho[n,m+1] + np.sqrt(n+1)*np.sqrt(m+1)*rho[n+1,m+1]
    
    def solve(self,n,m,rho_0):
        pass

    
if __name__ == '__main__':
    delta = 0
    epsilon = 1
    K = -0.5
    n_states = 3
    kr = Kerr_Resonator(n_states,delta,epsilon,K)
    rho_0 = (np.outer(kr.n_vector(0), kr.n_vector(0)))