import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as int

class Kerr_Resonator:


    def __init__(self, n_states):
        
        self.n_states = n_states    


    def lowering_operator(self):
        n = np.arange(1, self.n_states, 1)
        operator = np.diag(np.sqrt(n), k=1)
        return operator

    def raising_operator(self):
        n = np.arange(1, self.n_states, 1)
        operator = np.diag(np.sqrt(n), k=-1)
        return operator

    def n_vector(self, n):
        vector = np.zeros(self.n_states) 
        comp = np.arange(0,self.n_states,1)
        for i in comp:
            if i==n:
                vector[i-(self.n_states)] = 1
        return vector

    def rho_dot(self, t, rho):
        delta = 0
        epsilon = 0.5
        K = -0.5

        a = self.lowering_operator()
        a_dag = self.raising_operator()
        i = 1j #imaginary unit
        rho = rho.reshape(self.n_states, self.n_states)

        rho_d = - i*delta*(a@a_dag@a@rho) - i*(K/2)*(a@a_dag@a_dag@a@a@rho) - epsilon*(a@a_dag@rho)
        + epsilon*(a@a@rho)+ i*delta*(a@rho@a_dag@a) + i*(K/2)*(a@rho@a_dag@a_dag@a@a) - epsilon*(a@rho@a_dag)
        + epsilon*(a@rho@a) - (1/2)*(a_dag@a@rho - 2*(a@rho@a_dag) + rho@a_dag@a)

        return rho_d.flatten()
    
    def solve(self, rho_0, t_start, t_stop):

        soln = int.solve_ivp(self.rho_dot, [t_start, t_stop], rho_0.flatten())
        return soln


if __name__ == '__main__':
    kr = Kerr_Resonator(2)
    rho_0 = np.outer(kr.n_vector(0), kr.n_vector(0))
    soln = kr.solve(rho_0, 0, 1)

    print(soln.y)