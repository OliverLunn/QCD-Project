from qutip import *
import numpy as np
import matplotlib.pyplot as plt
N=10
rho_0 = fock_dm(N,0)
a = destroy(N)
a_dag = create(N)
print(a_dag)
K = -0.5
epsilon = 1
h_bar = 1.05e-34

deltas = np.arange(0,10,0.1)
photon_num = np.zeros(len(deltas))
g = np.zeros(len(deltas))
times = np.linspace(0, 20, 100)
i=0


for delta in deltas:
    H = (delta * a_dag * a + (K/2) * a_dag *a_dag * a * a + 1j * epsilon * (a_dag - a))
    result = mesolve(H, rho_0, times, [a], [a_dag*a, a_dag*a_dag*a*a])
    photon_num[i] = (result.expect[0])[-1]
    g[i] = (result.expect[1])[-1]/((result.expect[0])[-1])**2
    
    density_matrix = expect(result.states, projection(2,0,0))
    
    i+=1


plt.plot(deltas,photon_num,label='$<a^{+}a>$')
plt.plot(deltas,g,label='$g^{(2)} (0)$')
plt.legend()
plt.xlabel('$\Delta$')
plt.show()
print(density_matrix)


