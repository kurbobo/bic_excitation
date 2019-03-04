from matplotlib import pyplot as plt
import numpy as np
xi=1
Gamma=0.6# 6e-04
alpha= 0.8333338333333
delta_0=0.5
gamma_C=1
s=6
# B=np.arange(0.5,1,0.000001)
B=0.28430769429548375
g=2*xi*B*np.sqrt((gamma_C * delta_0 - gamma_C*alpha*np.abs(B)**2 + Gamma + s)**2 + (delta_0 - alpha*np.abs(B)**2 - Gamma * gamma_C )**2)
print(g)
# plt.plot(np.abs(g),B)
# plt.text(2,0.6,r"$\alpha=$"+str(alpha)+r"$ \gamma_C=$"+str(gamma_C)+r" s="+str(s))
#
# plt.xlabel('|g|')
# plt.ylabel('B_{20}')
# plt.show()
