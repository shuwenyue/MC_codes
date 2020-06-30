import os
import numpy as np
import matplotlib.pyplot as plt

file = 'Energies.11.txt'
density = np.genfromtxt(file, skip_header=0,invalid_raise=False)

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(density[:,0], density[:,1],linestyle='-',label='vapor')
ax1.plot(density[:,0], density[:,2],linestyle='-',label='liquid')

ax1.set_xlabel('g(r)')
ax1.set_ylabel('r (nm)')

fig.savefig('/home/syue/Projects/class/MSE/hw5/density.png', dpi=fig.dpi)

plt.show()


