from matplotlib import pyplot as plt
import numpy as np
def result():
    return np.sqrt(delta_light**2 - 3*gamma_light**2)
xi=1
s=16.48/16
delta=0.2
alpha=1
gamma=0.2
gamma_C=1
gamma_light=gamma + s/(1+gamma_C**2)
delta_light=delta + gamma_C*s/(1+gamma_C**2)
print('delta_light = ',delta_light)
print('gamma_light = ',gamma_light)
delta_0=np.arange(np.sqrt(3)*gamma,10*gamma,0.0001)
B_minus=np.sqrt((2*delta_0-np.sqrt(delta_0**2-3*gamma**2))/3/alpha)
B_plus=np.sqrt((2*delta_0+np.sqrt(delta_0**2-3*gamma**2))/3/alpha)
plt.plot(delta_0,B_minus, label="B_minus")
plt.plot(delta_0,B_plus, label="B_plus")
plt.ylabel('$|B_{20}|$')
plt.xlabel('$\delta_0$')
gamma_light=gamma + s/(1+gamma_C**2)
delta_light=delta_0 + gamma_C*s/(1+gamma_C**2)
B_plus_new= np.sqrt((2*delta_light +result())/3/alpha)
B_minus_new= np.sqrt((2*delta_light -result())/3/alpha)
plt.plot(delta_0,B_minus_new, label="gran_new")
plt.plot(delta_0,B_plus_new, label="gran2_new")
plt.grid(True)
plt.legend(loc='best')
plt.show()
# delta_0=0.7
# delta_light=delta_0 + gamma_C*s/(1+gamma_C**2)
# B_plus_new= np.sqrt((2*delta_light +result())/3/alpha)
# B_minus_new= np.sqrt((2*delta_light -result())/3/alpha)
# print(B_minus_new)
# print(B_plus_new)