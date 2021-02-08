import numpy as np
import common as co
import matplotlib.pyplot as plt 

h = 0.002

# planets
Venus = co.simplePlanetRK4(0.72,0.6152,h)
Earth = co.simplePlanetRK4(1,1,h)
Mars = co.simplePlanetRK4(1.52,1.8809,h)

# T^2
TV2 = (0.6152)**2
TE2 = (1)**2
TM2 = (1.8809)**2

T2 = np.array([TV2,TE2,TM2])

#a^3
aV3 = np.amax(co.dist(Venus[0],Venus[1]))**3
aE3 = np.amax(co.dist(Earth[0],Earth[1]))**3
aM3= np.amax(co.dist(Mars[0],Mars[1]))**3


a3 = np.array([aV3,aE3,aM3])

plt.plot(a3,T2,color='brown')
plt.plot(aV3,TV2, 'go', label='Venus')
plt.plot(aE3,TE2, 'bo', label='Earth')
plt.plot(aM3,TM2,'ro', label='Mars')
plt.xlabel('$a^3$', fontsize=14)
plt.ylabel('$T^2$', fontsize=14)
plt.legend(loc=2)
