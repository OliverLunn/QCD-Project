import numpy as np
import matplotlib.pyplot as plt
from qutip import *
from tqdm import tqdm

def matrix(n_states, result,t_steps):
    density_matrix = np.zeros((n_states, n_states,t_steps),dtype=complex)
    for n in range(N):
        for m in range(N):
            density_matrix[n,m,:] = (expect(result.states, projection(n_states, n, m)))
    return density_matrix

N=15
K= -0.5
K_prime = -0.5

times = np.linspace(0, 20, 200)
t_steps = len(times)

epsilons = np.linspace(1,5,10)
deltas = np.linspace(-15,15,10)

rho_0 = fock_dm(N,0)
a = destroy(N)
a_dag = create(N)

fig, ax = plt.subplots(1, figsize=(14, 7), facecolor='w', edgecolor='k')  #create subplot
axs = 0
P = np.zeros((len(deltas),len(epsilons)))
j=0
for epsilon in tqdm(epsilons):
    i=0
    for delta in tqdm(deltas):
        H = (delta * a_dag * a + (K/2) * a_dag *a_dag * a * a + 1j * epsilon * (a_dag - a) + (K_prime/3) * (a_dag**3) * (a**3))
        result = mesolve(H,rho_0,times,[a])
        density_mat = matrix(N,result,t_steps)
        P[i,j] = np.real(density_mat[1,1,-1])
        i+=1
    j+=1

ax.imshow(P,origin='lower',extent=(np.min(epsilons),np.max(epsilons),np.min(deltas),np.max(deltas)))
ax.set_ylabel('$\Delta_p$',fontsize="16")
ax.set_xlabel('$\epsilon_p$',fontsize="16")
#ax.legend()


plt.show()
plt.tight_layout()