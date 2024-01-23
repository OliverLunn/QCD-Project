from qutip import *
import numpy as np
import matplotlib.pyplot as plt

N=15
K = -5
K_prime = -5
epsilon = 1
times = np.linspace(0, 20, 550)
delta = 0

rho_0 = fock_dm(N,0)
a = destroy(N)
a_dag = create(N)

fig,ax = plt.subplots(1)
H = (delta * a_dag * a + (K/2) * a_dag *a_dag * a * a + 1j * epsilon * (a_dag - a) + (K_prime/3) * (a_dag**3) * (a**3))
result = mesolve(H,rho_0,times,[a],[a_dag*a,a_dag*a_dag*a*a])
photon_num = result.expect[0]

ax.plot(times,result.expect[0],label='$K=$'+str(K)+', $K^{\prime}=$'+str(K_prime))
ax.set_xlabel('t')
ax.set_ylabel('Photon number $<a^+ a>$')

plt.legend()
plt.show()