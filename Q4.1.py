import numpy as np
import common as co
import matplotlib.pyplot as plt 


EarthMars = co.threeBodySystemRK4(1,1,1.52,1.8809,0.001,10)
#EarthMars = co.threeBodySystemRK4(1,1,1,2,2.2,0.1074,0.001,8)

plt.plot(EarthMars[0],EarthMars[1],linewidth = 0.5, label='Earth numerical orbit',color='navy')
plt.plot(EarthMars[2],EarthMars[3],linewidth = 1, label='Mars numerical orbit',color='crimson')
plt.plot(0,0, 'yo',markersize=15,label='Sun', color='gold')
plt.xlabel(r"$x \,[\mathrm{AU}]$", fontsize=14)
plt.ylabel(r"$y \,[\mathrm{AU}]$", fontsize=14)
plt.grid()
plt.legend(loc=1)
