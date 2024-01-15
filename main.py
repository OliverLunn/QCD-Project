import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import odeintw

class Kerr_Resonator:

    def __init__(self, n_states):
        
        self.n_states = n_states    

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
        delta = 0
        epsilon = 0.5
        K = -0.5

        a = self.lowering_operator()
        a_dag = self.raising_operator()
        i = 1j #imaginary unit

        rho_d = - i*delta*np.dot(a@a_dag@a,rho) - i*(K/2)*np.dot(a@a_dag@a_dag@a@a,rho) - epsilon*np.dot(a@a_dag,rho)
        + epsilon*np.dot(a@a,rho)+ i*delta*np.dot(np.dot(a,rho),a_dag@a) + i*(K/2)*np.dot(np.dot(a,rho),a_dag@a_dag@a@a) - epsilon*np.dot(np.dot(a,rho),a_dag)
        + epsilon*np.dot(np.dot(a,rho),a) - (1/2)*(np.dot(a_dag@a,rho) - 2*np.dot(a,rho)@a_dag) + np.dot(rho,a_dag@a)
    
        return rho_d
    
    def solve(self, rho_0, t_start, t_stop):

        soln = odeintw.odeintw(self.rho_dot,rho_0,np.linspace(t_start,t_stop,1000))
        return soln


if __name__ == '__main__':
    kr = Kerr_Resonator(3)
    t=0
    rho_0 = (np.outer(kr.n_vector(0), kr.n_vector(0)))
    soln = kr.solve(rho_0, 0, 100)
    plt.plot(np.real(soln[:,1,1]),'.',color = 'blue')
    plt.plot(np.real(soln[:,0,0]),'.',color = 'red')
    plt.plot(np.real(soln[:,0,1]),'.',color = 'green')
    plt.plot(np.real(soln[:,1,0]),'.',color = 'orange')
    plt.show()