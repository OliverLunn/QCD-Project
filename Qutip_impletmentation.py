from qutip import *
import numpy as np
import matplotlib.pyplot as plt
N=10
rho_0 = fock_dm(N,0)
a = destroy(N)
a_dag = create(N)

K = 0
epsilon = 1


times = np.linspace(0, 20, 100)
delta = 0

fig,ax = plt.subplots(1)
H = (delta * a_dag * a + (K/2) * a_dag *a_dag * a * a + 1j * epsilon * (a_dag - a))
result = mesolve(H,rho_0,times,[a],[a_dag*a,a_dag*a_dag*a*a])

ax.plot(times,result.expect[0],label='$<a^{+}a>$')
ax.set_xlabel('t')
ax.set_ylabel('expctation value')

fig,ax2 = plt.subplots(1)
deltas = np.arange(-10,10,0.1)
photon_num = np.zeros(len(deltas))
g = np.zeros(len(deltas))

i=0
for delta in deltas:
    H = (delta * a_dag * a + (K/2) * a_dag *a_dag * a * a + 1j * epsilon * (a_dag - a))
    result = mesolve(H,rho_0,times,[a],[a_dag*a,a_dag*a_dag*a*a])
    photon_num[i] = (result.expect[0])[-1]
    g[i] = (result.expect[1])[-1]/((result.expect[0])[-1])**2
    i+=1


ax2.plot(deltas,photon_num,label='$<a^{+}a>$')

#ax2.plot(deltas,g,label='$g^{(2)} (0)$')
ax2.legend()
ax2.set_xlabel('$\Delta$')

plt.show()