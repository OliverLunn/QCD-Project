import numpy as np
import matplotlib.pyplot as plt
from qutip import *
from tqdm import tqdm

N=10
rho_0 = fock_dm(N,0)
a = destroy(N)
a_dag = create(N)

epsilon = 1
delta = 1
K= 1
K_prime = 0

times = np.linspace(0, 20, 1000)

H = (delta * a_dag * a + (K/2) * a_dag *a_dag * a * a + 1j * epsilon * (a_dag - a) + (K_prime/3) * (a_dag**3) * (a**3))
result = mesolve(H,rho_0,times,[a],[a_dag*a,a_dag*a_dag*a*a])

plt.plot(times,result.expect[0])
plt.show()