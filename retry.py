import numpy as np
import matplotlib.pyplot as plt
import odeintw

def rho_dot_element(n:int,m:int,rho,delta,epsilon,K):
    i=1j
    rho = np.pad(rho,1,'constant',constant_values=0+0j)
    rho_dot = (-i*delta*n - i*(K/2)*n*(n-1) + i*delta*m + i*(K/2)*m*(m-1) - (1/2)*n - (1/2)*m)*rho[n,m] 
    + epsilon*np.sqrt(n)*rho[n-1,m] + epsilon*np.sqrt(m)*rho[n,m-1] - epsilon*np.sqrt(n+1)*rho[n+1,m] 
    - epsilon*np.sqrt(m+1)*rho[n,m+1] + np.sqrt(n+1)*np.sqrt(m+1)*rho[n+1,m+1]
    return rho_dot

def n_vec(n_states:int,n:int):
    a = np.arange(0,n_states,1)
    vec = np.where(a==n,1+0j,0+0j)
    return vec

def rho_dot_matrix(rho,t,n_states:int,delta,epsilon,K):
    rho_dot = np.empty((n_states,n_states),dtype=complex)
    for i in range(n_states):
        for j in range(n_states):
            rho_dot[i,j] = rho_dot_element(i,j,rho,delta,epsilon,K)
    return rho_dot

def photon_number(rho, t_idx,n_states):
    """
    Calculates the number of photons in the oscillator. 
    """
    a = raising_operator(n_states)
    a_dag = lowering_operator(n_states)
    expectation_val = np.trace(a_dag@a@rho[t_idx,:,:])

    return expectation_val

def lowering_operator(n_states):
    n = np.arange(1, n_states, 1)
    operator = np.diag(np.sqrt(n), k=1)
    return operator

def raising_operator(n_states):
    n = np.arange(1, n_states, 1)
    operator = np.diag(np.sqrt(n), k=-1)
    return operator

n_states = 5
times = np.linspace(0, 20, 100)
delta = 0
epsilon = 1
K = -0.5
rho_0 = np.outer(n_vec(n_states,0),n_vec(n_states,0))
result = odeintw.odeintw(rho_dot_matrix,rho_0,times,args=(n_states,delta,epsilon,K))


p_num = np.empty(len(times))
for i in range(len(times)):
    p_num[i] = photon_number(result,i,n_states)
plt.plot(times,p_num)
plt.show()