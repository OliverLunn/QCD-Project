from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

def alpha_sq(epsilon,delta):
    return ((epsilon)**2 )/((delta)**2 + 1/4)

N=15
K = 0
epsilon = 1
times = np.linspace(0, 20, 100)

rho_0 = fock_dm(N,0)
a = destroy(N)
a_dag = create(N)

fig,ax1 = plt.subplots(1)
deltas = np.arange(-5,5,0.05)
photon_num = np.zeros(len(deltas))
photon_num_analytic = np.zeros(len(deltas))
g = np.zeros(len(deltas))

i=0
for delta in tqdm(deltas):
    H = (delta * a_dag * a + (K/2) * a_dag *a_dag * a * a + 1j * epsilon * (a_dag - a))
    result = mesolve(H,rho_0,times,[a],[a_dag*a,a_dag*a_dag*a*a])
    photon_num[i] = (result.expect[0])[-1]
    photon_num_analytic[i] = alpha_sq(epsilon,delta)
    g[i] = (result.expect[1])[-1]/((result.expect[0])[-1])**2
    i+=1

ax1.plot(deltas, photon_num, ".-b", label='numeric')
ax1.plot(deltas, photon_num_analytic, "D-k", label='analytic')
ax1.set_xlabel('$\Delta$',fontsize="16")
ax1.set_ylabel('$<a^{+}a>$',fontsize="16")
ax1.legend()
plt.show()