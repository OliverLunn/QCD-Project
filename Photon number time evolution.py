from qutip import *
import numpy as np
import matplotlib.pyplot as plt

N=10
K = -0.5
epsilon = 1
times = np.linspace(0, 20, 100)
delta = 0

rho_0 = fock_dm(N,0)
a = destroy(N)
a_dag = create(N)

fig,ax = plt.subplots(1)
H = (delta * a_dag * a + (K/2) * a_dag *a_dag * a * a + 1j * epsilon * (a_dag - a))
result = mesolve(H,rho_0,times,[a],[a_dag*a,a_dag*a_dag*a*a])

ax.plot(times,result.expect[0],label='$<a^{+}a>$')
ax.set_xlabel('t')
ax.set_ylabel('Photon number $<a^+ a>$')
plt.show()
