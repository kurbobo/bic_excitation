from matplotlib import pyplot as plt
import numpy as np
def result():
    return np.sqrt(delta_light**2 - 3*gamma_light**2)
xi=1
delta=0.6
alpha=1
gamma=0.1
gamma_C=1e-1
s=3e-01
B_minus=np.sqrt((2*delta-np.sqrt(delta**2-3*gamma**2))/3/alpha)
B_plus=np.sqrt((2*delta+np.sqrt(delta**2-3*gamma**2))/3/alpha)
print(B_minus, " ",B_plus)
delta_0=np.arange(np.sqrt(3)*gamma,10*gamma,0.0001)
B_minus=np.sqrt((2*delta_0-np.sqrt(delta_0**2-3*gamma**2))/3/alpha)
B_plus=np.sqrt((2*delta_0+np.sqrt(delta_0**2-3*gamma**2))/3/alpha)
plt.plot(delta_0,B_minus, label="B_minus")
plt.plot(delta_0,B_plus, label="B_plus")
a=gamma_C*delta+gamma+s
b=delta-gamma_C*gamma
c=gamma_C*alpha
res=alpha**2*(gamma_C**2+1)**2*(delta_0**2-3*gamma**2)+2*s*(gamma_C**2+1)*alpha**2* (gamma_C*delta_0-3*gamma)+s**2*alpha**2*(gamma_C**2-3)
B_pl= np.sqrt((np.sqrt(res)+2*(a*c+ alpha * b))/3/(alpha**2 + c**2))
B_min= np.sqrt((-np.sqrt(res)+2*(a*c+ alpha * b))/3/(alpha**2 + c**2))
plt.plot(delta_0,B_min, label="gran")
plt.plot(delta_0,B_pl, label="gran2")
plt.ylabel('$|B_{20}|$')
plt.xlabel('$\delta_0$')
a=gamma_C*delta_0+gamma+s
b=delta_0-gamma_C*gamma
gamma_light=gamma + s/(1+gamma_C**2)
delta_light=delta_0 + gamma_C*s/(1+gamma_C**2)
B_plus_new= np.sqrt((2*delta_light +result())/3/alpha)
B_minus_new= np.sqrt((2*delta_light -result())/3/alpha)
plt.plot(delta_0,B_minus_new, label="B_minus_new")
plt.plot(delta_0,B_plus_new, label="B_plus_new")
c=gamma_C*alpha
B_plus= np.sqrt(np.sqrt(4*( a*c + alpha *b)**2 - (a**2 + b**2) *(3 * alpha**2 + 3*c**2))+2*(a*c+ alpha * b)/3/(alpha**2 + c**2))
B_minus= np.sqrt(-np.sqrt(4*( a*c + alpha *b)**2 - (a**2 + b**2) *(3 * alpha**2 + 3*c**2))+2*(a*c+ alpha * b)/3/(alpha**2 + c**2))
plt.grid(True)
plt.legend(loc='best')
plt.show()
delta_0=0.7
delta_light=delta_0 + gamma_C*s/(1+gamma_C**2)
B_plus_new= np.sqrt((2*delta_light +result())/3/alpha)
B_minus_new= np.sqrt((2*delta_light -result())/3/alpha)
print(B_minus_new)
print(B_plus_new)