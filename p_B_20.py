from matplotlib import pyplot as plt
import numpy as np
def g(B):
    return 2*xi*B*np.sqrt((gamma_C * delta - gamma_C*alpha*np.abs(B)**2 + gamma + s)**2 + (delta - alpha*np.abs(B)**2 - gamma * gamma_C )**2)
def re(gamma_C):
    return alpha**2*(gamma_C**2+1)**2*(delta**2-3*gamma**2)+2*s*(gamma_C**2+1)*alpha**2* (gamma_C*delta-3*gamma)+s**2*alpha**2*(gamma_C**2-3)
# def real_func(B_list):
#     g_r=[]
#     for B in B_list:
#         if B>=0.196269 and B<=0.712119: g_r.append(14.4276)
#         else: g_r.append(g(B))
#     return g_r
xi=1
delta=0.2#0.6
alpha=1
gamma=6e-02
gamma_C=1e-01
s=1.648e+01
B_min=np.sqrt(delta/3/alpha)
B_max=np.sqrt(delta/alpha)
print("g_min=",g(B_min))
print("g_max=",g(B_max))
B=np.arange(0,0.76,0.000001)
gamma_Ce=np.arange(0,0.76,0.000001)
plt.plot(B,np.abs(g(B)))
# for i in B:
#     if np.abs(g(i)- 14.4276)<0.000035:
#         print(i)
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
plt.ylabel('|g|')
plt.xlabel('B_{20}')
plt.show()
