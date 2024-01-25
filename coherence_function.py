from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

def alpha_sq(epsilon,delta):
    return ((epsilon)**2 )/((delta)**2 + 1/4)


def g(deltas, K, K_prime):

    g = np.zeros(len(deltas))
    i=0
    for delta in tqdm(deltas):
        H = (delta * a_dag * a + (K/2) * a_dag *a_dag * a * a + 1j * epsilon * (a_dag - a) + (K_prime/3) * (a_dag**3) * (a**3))
        result = mesolve(H,rho_0,times,[a],[a_dag*a, a_dag*a_dag*a*a])
        g[i] = (result.expect[1])[-1]/(((result.expect[0])[-1])*np.conj(result.expect[0])[-1])
        i+=1
    return np.real(g)

N = 15
epsilon = 1
times = np.linspace(0, 20, 100)
deltas = np.arange(-8,20,0.1)

rho_0 = fock_dm(N,0)
a = destroy(N)
a_dag = create(N)

fig, ax = plt.subplots(1,1, figsize=(14, 7), facecolor='w', edgecolor='k')  #create subplot


ax.plot(deltas, g(deltas, 0, 0), ".-k", label='$^{\'}K/\kappa$=0')
ax.plot(deltas, g(deltas, 0, -1), ".-r", label='$K^{\'}/\kappa$=-1')
ax.plot(deltas, g(deltas, 0, -2), ".-g", label='$K^{\'}/\kappa$=-2')

ax.set_xlabel('$\Delta_p/\kappa$', fontsize="32")
ax.set_ylabel('$g^2(0)$', fontsize="32")
ax.tick_params(axis="x", labelsize=28)
ax.tick_params(axis="y", labelsize=28)
plt.legend(fontsize="24")
plt.tight_layout()
plt.show()