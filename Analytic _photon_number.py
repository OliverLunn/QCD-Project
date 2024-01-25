from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from qutip import *
from tqdm import tqdm
import numpy.ma as ma

epsilon = 1
K= -0.5

alpha = np.arange(0,1000,0.1)
times = np.linspace(0, 20, 100)

N=10
rho_0 = fock_dm(N,0)
a = destroy(N)
a_dag = create(N)

fig,ax1 = plt.subplots(1)
deltas = np.arange(-5,5,0.05)
photon_num = np.zeros(len(deltas))
roots = np.zeros((len(deltas),3),dtype=complex)

i=0
for delta in tqdm(deltas):
    poly= np.polynomial.polynomial.Polynomial([-epsilon**2,(delta**2 + 1/4),2*delta*K,K**2])
    roots[i] = poly.roots()

    H = (delta * a_dag * a + (K/2) * a_dag *a_dag * a * a + 1j * epsilon * (a_dag - a))
    result = mesolve(H,rho_0,times,[a],[a_dag*a,a_dag*a_dag*a*a])
    photon_num[i] = (result.expect[0])[-1]
    i+=1

roots = ma.masked_where(np.isreal(roots)==False,roots)

#photon_num_analytic = [i for i in photon_num_analytic if i.imag==0]
for i in range(3):
    ax1.plot(deltas,roots[:,i],'x',label='Semi Classical Approximation')
ax1.plot(deltas, photon_num, "x-r", label='numeric')
ax1.set_xlabel('$\Delta$',fontsize="16")
ax1.set_ylabel('$<a^{+}a>$',fontsize="16")
ax1.legend()
plt.show()
