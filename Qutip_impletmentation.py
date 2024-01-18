from qutip import *
import numpy as np
import matplotlib.pyplot as plt

def density_matrix(n_states):
    density_matrix = np.zeros((n_states, n_states))

    for i in range(n_states):
        for j in range(n_states):
            density_matrix[i,j] = expect(result.states, projection(n_states, i, j))

    return density_matrix

N=3
rho_0 = fock_dm(N,0)
a = destroy(N)
a_dag = create(N)
K = -0.5
epsilon = 1
h_bar = 1.05e-34

deltas = np.arange(-5,5,0.01)
photon_num = np.zeros(len(deltas))
g = np.zeros(len(deltas))
times = np.linspace(0, 20, 250)
i=0


for delta in deltas:
    H = (delta * a_dag * a + (K/2) * a_dag *a_dag * a * a + 1j * epsilon * (a_dag - a))
    result = mesolve(H, rho_0, times, [a], [a_dag*a, a_dag*a_dag*a*a])
    photon_num[i] = (result.expect[0])[-1]
    g[i] = (result.expect[1])[-1]/((result.expect[0])[-1])**2
    
    density_matrix = expect(result.states, projection(2,0,0))
    
    i+=1

#density_mat = density_matrix(N)
plt.plot(deltas,photon_num,label='$<a^{+}a>$')
plt.plot(deltas,g,label='$g^{(2)} (0)$')
plt.legend()
plt.xlabel('$\Delta$')
plt.show()
print(density_matrix)


