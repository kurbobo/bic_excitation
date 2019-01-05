from matplotlib import pyplot as plt
import numpy as np
Gamma=0.01
sigma=0.01
sigma2=0.1
# delta_0=1
alpha=0.01
delta_0=np.arange(np.sqrt(3/4)*Gamma,15*Gamma,0.001)
B_minus=np.sqrt((4*delta_0-2*np.sqrt(4*delta_0**2-3*Gamma**2))/alpha)
B_plus=np.sqrt((4*delta_0+2*np.sqrt(4*delta_0**2-3*Gamma**2))/alpha)

plt.plot(delta_0,B_minus)
plt.plot(delta_0,B_plus)
plt.ylabel('$|B_{20}|$')
plt.xlabel('$\delta_0$')
plt.show()
delta_0=3
B=np.sqrt((4*delta_0+2*np.sqrt(4*delta_0**2-3*Gamma**2))/alpha)
print(B)
