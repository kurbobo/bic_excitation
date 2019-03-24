from matplotlib import pyplot as plt
import numpy as np
def p(B):
    return  B*np.sqrt(gamma_light**2+(delta_light - alpha*B**2)**2)
def result():
    return np.sqrt(delta_light**2 - 3*gamma_light**2)
xi=1
s=0.1
delta=0.6
alpha=1
gamma=0.1
gamma_C=0.1
gamma_light=gamma + s/(1+gamma_C**2)
delta_light=delta + gamma_C*s/(1+gamma_C**2)
dist=1000
B=np.arange(0,0.9,0.000001)
plt.plot(B,np.abs(p(B)))
# res=alpha**2*(gamma_C**2+1)**2*(delta**2-3*gamma**2)+2*s*(gamma_C**2+1)*alpha**2* (gamma_C*delta-3*gamma)+s**2*alpha**2*(gamma_C**2-3)
B_plus= np.sqrt((2*delta_light +result())/3/alpha)
B_minus= np.sqrt((2*delta_light -result())/3/alpha)
plt.axvline(x=B_plus,color='green')
plt.axvline(x=B_minus,color='red')
plt.ylabel('|p|')
plt.xlabel('B_{20}')
plt.show()
