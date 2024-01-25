from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

def matrix(n_states, result,t_steps):

    density_matrix = np.zeros((n_states, n_states,t_steps),dtype=complex)
    for n in range(N):
        for m in range(N):
            density_matrix[n,m,:] = (expect(result.states, projection(n_states, n, m)))
    return density_matrix
        
            
N = 15
K = 0
K_prime = 0
epsilon = 1
times = np.linspace(0, 20, 200)
t_steps = len(times)
deltas = np.linspace(-10,10,201)

trace = np.zeros(len(deltas))

rho_0 = fock_dm(N,0)
a = destroy(N)
a_dag = create(N)

fig = plt.figure()
ax = plt.axes([0.1,0.15,0.8,0.8])
ax1 = plt.axes([0.15,0.2,0.25,0.25])
i=0
for delta in tqdm(deltas):
    H = (delta * a_dag * a + (K/2) * a_dag *a_dag * a * a + 1j * epsilon * (a_dag - a) + (K_prime/3)*a_dag*a_dag*a_dag*a*a*a)
    result = mesolve(H, rho_0, times,[a])
    density_m = matrix(N, result,t_steps)
    trace[i] = np.real(np.trace(density_m[:,:,-1]@np.conj(density_m[:,:,-1])))
    i+=1

ax.plot(deltas, trace, ".-b")
ax.set_ylabel("$Tr[\\rho^2]$",fontsize="35")
ax.set_xlabel("$\Delta_p / \kappa $", fontsize="35")

deltas2 = np.linspace(-1.0,1.0,201)
trace2 = np.zeros(len(deltas2))
i=0
for delta in tqdm(deltas2):
    H = (delta * a_dag * a + (K/2) * a_dag *a_dag * a * a + 1j * epsilon * (a_dag - a) + (K_prime/3)*a_dag*a_dag*a_dag*a*a*a)
    result = mesolve(H, rho_0, times,[a])
    density_m = matrix(N, result,t_steps)
    trace2[i] = np.real(np.trace(density_m[:,:,-1]@np.conj(density_m[:,:,-1])))
    i+=1

ax1.plot(deltas2, trace2, ".-b")

ax.tick_params('both',labelsize='24')
plt.show()


