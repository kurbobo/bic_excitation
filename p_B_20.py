from matplotlib import pyplot as plt
import numpy as np
def g(B):
    return 2*xi*B*np.sqrt((gamma_C * delta - gamma_C*alpha*np.abs(B)**2 + gamma + s)**2 + (delta - alpha*np.abs(B)**2 - gamma * gamma_C )**2)
def re(gamma_C):
    return alpha**2*(gamma_C**2+1)**2*(delta**2-3*gamma**2)+2*s*(gamma_C**2+1)*alpha**2* (gamma_C*delta-3*gamma)+s**2*alpha**2*(gamma_C**2-3)
xi=1
s=0.1
delta=0.6
alpha=1
gamma=0.1
gamma_C=0.1
B_min=np.sqrt(delta/3/alpha)
B_max=np.sqrt(delta/alpha)
print("g_min=",str(g(B_min)) + " B_min = " + str(B_min))
print("g_max=",str(g(B_max))+ " B_max = " + str(B_max))
B=np.arange(0,0.9,0.000001)
# gamma_Ce=np.arange(0,0.76,0.000001)
plt.plot(B,np.abs(1j*g(B)/2/xi/(1j+gamma_C)))
plt.text(0.0,23,r"$\alpha=$ "+str(alpha)+r" $ \gamma_C=$"+str(gamma_C)+r" s="+str(s)+r" $ \delta=$"+str(delta))
a=gamma_C*delta+gamma+s
b=delta-gamma_C*gamma
c=gamma_C*alpha
res=alpha**2*(gamma_C**2+1)**2*(delta**2-3*gamma**2)+2*s*(gamma_C**2+1)*alpha**2* (gamma_C*delta-3*gamma)+s**2*alpha**2*(gamma_C**2-3)
B_plus= np.sqrt((np.sqrt(res)+2*(a*c+ alpha * b))/3/(alpha**2 + c**2))
B_minus= np.sqrt((-np.sqrt(res)+2*(a*c+ alpha * b))/3/(alpha**2 + c**2))
plt.axvline(x=B_plus,color='green')
plt.axvline(x=B_minus,color='red')
# plt.plot(gamma_Ce,re(gamma_Ce))
plt.ylabel('|p|')
plt.xlabel('B_{20}')
plt.show()
gamma_light=gamma + s/(1+gamma_C**2)
delta_light=delta + gamma_C*s/(1+gamma_C**2)
print('delta_light = ',delta_light)
print('gamma_light = ',gamma_light)
p=1j*g/2/xi/(1j+gamma_C)
print('p = ',p)
