import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from numpy import vectorize


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

    def rho_dot(self, t, rho, row):
        delta = 0
        epsilon = 0.5
        K = -0.5

        a = self.lowering_operator()
        a_dag = self.raising_operator()
        i = 1j #imaginary unit

        rho_d = - i*delta*(a@a_dag@a@rho) - i*(K/2)*(a@a_dag@a_dag@a@a@rho) - epsilon*(a@a_dag@rho)
        + epsilon*(a@a@rho)+ i*delta*(a@rho@a_dag@a) + i*(K/2)*(a@rho@a_dag@a_dag@a@a) - epsilon*(a@rho@a_dag)
        + epsilon*(a@rho@a) - (1/2)*(a_dag@a@rho - 2*(a@rho@a_dag) + rho@a_dag@a)

        return rho_d[row]
    
    def solve(self, rho_0, t_start, t_stop):
        rows = self.n_states
        t_points = np.linspace(t_start, t_stop, 4)
        soln_matrix = np.zeros((rows, rows, len(t_points)))

        for r in range(rows):
                soln = integrate.solve_ivp(self.rho_dot, [t_start, t_stop], rho_0[r,:], t_eval=t_points, args=(r,))
                #print(soln.y.shape)
                soln_matrix[r,:,:] = soln.y
                print(soln_matrix)

        return soln_matrix


if __name__ == '__main__':
    kr = Kerr_Resonator(3)
    t_start, t_stop = 0, 1
    rho_0 = np.outer(kr.n_vector(0), kr.n_vector(0))
    soln = kr.solve(rho_0, t_start, t_stop)