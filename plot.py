import numpy as np
import matplotlib.pyplot as plt
P = np.loadtxt('P(1).dat')
fig = plt.figure()
ax = plt.axes([0.1,0.15,0.7,0.7])
epsilons = np.linspace(1,5,100)
deltas = np.linspace(-15,15,100)

plt.imshow(P,origin='lower',aspect='auto',extent=(np.min(epsilons),np.max(epsilons),np.min(deltas),np.max(deltas)))
ax.set_ylabel('$\Delta_p$',fontsize="16")
ax.set_xlabel('$\epsilon_p$',fontsize="16")

cax = plt.axes([0.85,0.15,0.05,0.7])
plt.colorbar(cax=cax)
cax.set_ylabel('$P(|1>)$')


plt.show()
np.savetxt('P(1).dat',P,fmt='%.4e')