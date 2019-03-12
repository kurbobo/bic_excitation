from matplotlib import pyplot as plt
import numpy as np
delta=0.2#0.6
alpha=1
gamma=6e-02
gamma_C=1e-01
s=6e+01
delta_0=0.2
B_minus=np.sqrt((2*delta_0-np.sqrt(delta_0**2-3*gamma**2))/3/alpha)
B_plus=np.sqrt((2*delta_0+np.sqrt(delta_0**2-3*gamma**2))/3/alpha)
print(B_minus, " ",B_plus)
delta_0=np.arange(np.sqrt(3)*gamma,10*gamma,0.001)
B_minus=np.sqrt((2*delta_0-np.sqrt(delta_0**2-3*gamma**2))/3/alpha)
B_plus=np.sqrt((2*delta_0+np.sqrt(delta_0**2-3*gamma**2))/3/alpha)
plt.plot(delta_0,B_minus, label="B_minus")
plt.plot(delta_0,B_plus, label="B_plus")

plt.ylabel('$|B_{20}|$')
plt.xlabel('$\delta_0$')
a=gamma_C*delta_0+gamma+s
b=delta_0-gamma_C*gamma
c=gamma_C*alpha
B_plus= np.sqrt(np.sqrt(4*( a*c + alpha *b)**2 - (a**2 + b**2) *(3 * alpha**2 + 3*c**2))+2*(a*c+ alpha * b)/3/(alpha**2 + c**2))
B_minus= np.sqrt(-np.sqrt(4*( a*c + alpha *b)**2 - (a**2 + b**2) *(3 * alpha**2 + 3*c**2))+2*(a*c+ alpha * b)/3/(alpha**2 + c**2))
plt.grid(True)
plt.legend(loc='best')
plt.show()
# delta_0=1
# B=np.sqrt((4*delta_0+2*np.sqrt(4*delta_0**2-3*gamma**2))/alpha)
# print(B)