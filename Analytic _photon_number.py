from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from qutip import *
from tqdm import tqdm

def interpolated_intercept(x, y1, y2):
    """Find the intercept of two curves, given by the same x data"""

    def intercept(point1, point2, point3, point4):
        """find the intersection between two lines
        the first line is defined by the line between point1 and point2
        the first line is defined by the line between point3 and point4
        each point is an (x,y) tuple.

        So, for example, you can find the intersection between
        intercept((0,0), (1,1), (0,1), (1,0)) = (0.5, 0.5)

        Returns: the intercept, in (x,y) format
        """    

        def line(p1, p2):
            A = (p1[1] - p2[1])
            B = (p2[0] - p1[0])
            C = (p1[0]*p2[1] - p2[0]*p1[1])
            return A, B, -C

        def intersection(L1, L2):
            D  = L1[0] * L2[1] - L1[1] * L2[0]
            Dx = L1[2] * L2[1] - L1[1] * L2[2]
            Dy = L1[0] * L2[2] - L1[2] * L2[0]

            x = Dx / D
            y = Dy / D
            return x,y

        L1 = line([point1[0],point1[1]], [point2[0],point2[1]])
        L2 = line([point3[0],point3[1]], [point4[0],point4[1]])

        R = intersection(L1, L2)

        return R

    idx = np.argwhere(np.diff(np.sign(y1 - y2)) != 0)
    xc, yc = intercept((x[idx], y1[idx]),((x[idx+1], y1[idx+1])), ((x[idx], y2[idx])), ((x[idx+1], y2[idx+1])))
    return xc,yc


def func1(alpha):
    return alpha*np.conj(alpha)

def func2(alpha,epsilon,K,delta):
    return (epsilon**2) / ((delta + K* (alpha*np.conj(alpha)))**2 + 1/4)

epsilon = 1
delta = 1
K= -0.5

alpha = np.arange(0,10,0.01)
f = func1(alpha)
g = func2(alpha,epsilon,K,delta)

alpha_c, fc = interpolated_intercept(alpha,f,g)
fc = fc[0][0]

times = np.linspace(0, 20, 100)

N=10
rho_0 = fock_dm(N,0)
a = destroy(N)
a_dag = create(N)

fig,ax1 = plt.subplots(1)
deltas = np.arange(-5,5,0.05)
photon_num = np.zeros(len(deltas))
photon_num_analytic = np.zeros(len(deltas))


i=0
for delta in tqdm(deltas):
    photon_num_analytic[i] = func2(float(fc),epsilon,K,delta)
    H = (delta * a_dag * a + (K/2) * a_dag *a_dag * a * a + 1j * epsilon * (a_dag - a))
    result = mesolve(H,rho_0,times,[a],[a_dag*a,a_dag*a_dag*a*a])
    photon_num[i] = (result.expect[0])[-1]
    g[i] = (result.expect[1])[-1]/((result.expect[0])[-1])**2
    i+=1

ax1.plot(deltas,photon_num_analytic,".-r", label='semi-clasical')
ax1.plot(deltas, photon_num, ".-b", label='numeric')
ax1.set_xlabel('$\Delta$',fontsize="16")
ax1.set_ylabel('$<a^{+}a>$',fontsize="16")
ax1.legend()
plt.show()
