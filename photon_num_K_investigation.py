from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

def alpha_sq(epsilon,delta):
    return ((epsilon)**2 )/((delta)**2 + 1/4)

def photon(deltas, K):
    photon_n = np.zeros(len(deltas))
    i=0
    for delta in tqdm(deltas):
        H = (delta * a_dag * a + (K/2) * a_dag *a_dag * a * a + 1j * epsilon * (a_dag - a))
        result = mesolve(H,rho_0,times,[a],[a_dag*a,a_dag*a_dag*a*a])
        photon_n[i] = result.expect[0][-1]
        i+=1
    return photon_n

N=15
epsilon = 1
times = np.linspace(0, 20, 100)
deltas = np.arange(-10,10, 0.1)

rho_0 = fock_dm(N,0)
a = destroy(N)
a_dag = create(N)

fig, ax = plt.subplots(1,1, figsize=(14, 7), facecolor='w', edgecolor='k')  #create subplot


ax.plot(deltas, photon(deltas, 0), ".-b", label='$K=0$')
ax.plot(deltas, photon(deltas, -1), ".-g", label='$K=-1$')
ax.plot(deltas, photon(deltas, -5), ".-r", label='$K=-5$')
ax.plot(deltas, photon(deltas, -10), ".-k", label='$K=-10$')
ax.set_xlabel('$\Delta$',fontsize="16")
ax.set_ylabel('$<a^{+}a>$',fontsize="16")

plt.legend()
plt.show()
plt.tight_layout()