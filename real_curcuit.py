import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
finite_time=50#should be less then 50
delta_0=1.1
chi_21=5
chi_31=17.11111
chi_22=chi_21
chi_32=chi_31
R=50/2.04
p=3.145
r_1=50
r_2=r_1
R_1=1
R_2=R_1
C_1=2E-7
C_2=C_1
C=1.1E-7
L_1=0.5E-3
L_2=L_1
w=np.sqrt(L_1/C_1)
k=1/np.sqrt(L_1*C_1)
sigma=w/R
varepsilon_0=p/sigma/np.sqrt(2)
print('w=',w)
print('k=',k)

# def alpha(omega):
#     return omega*(omega**2*(12*chi_31+2*chi_21**2)-3*chi_31**2)/(4*(4*omega**2-1))
alpha=3*chi_31/4+chi_21**2/6
print('alpha=',alpha)
# t=range(1,100000,1)
# p=[]
# for i in range(1,100000,1):p.append(alpha(i))
# plt.plot(t,p,'g',label='alpha')
# plt.show()
# lim=(w**2+R_1*r_1)/(2*r_1*w)/np.sqrt(1-4*alpha)
# print('B_20^2 must be >',lim)
# function that returns dz/dt
def model(z,t,varepsilon):
    noise =0# random.random()/10.0
    U_1 = z[0]+noise
    V_1 = z[1]
    U_2 = z[2]+noise
    V_2 = z[3]
    U_C = z[4]
    V_C = w / R*(U_1 - U_2-U_C- varepsilon)
    dU_1dt = (1 + chi_21 * U_1 + chi_31 * U_1** 2)*(-w / r_1 * U_1 + V_1 - V_C)
    dV_1dt = (-R_1 / w * z[1] - z[0])
    dU_2dt = (1 + chi_22 * z[2] + chi_32 * z[2] ** 2)*(-w * C_1 / r_2 / C_2 * z[2] + C_1 / C_2 * V_C - C_1 / C_2 * z[3])
    dV_2dt = (-R_2 * L_1 / w / L_2 * z[3] + L_1 / L_2 * z[2])
    dU_Cdt =w*C_1*(-1 / R / C * z[4] + 1 / C / R * (z[0] - z[2]) - 1 / C / R * varepsilon)
    dzdt = [dU_1dt, dV_1dt, dU_2dt, dV_2dt, dU_Cdt]
    return dzdt

# initial condition
z_0 = [0.01,0,0.0,0,0]

# number of time points
n = 500

# time points
dist=50
t = np.linspace(0,dist,n)

# step input
varepsilon = varepsilon_0*np.cos((1+np.cos(delta_0*t))*t)*np.heaviside(20-t,0)

# store solution
U_1 =np.empty_like(t)
V_1 = np.empty_like(t)
U_2 = np.empty_like(t)
V_2 = np.empty_like(t)
U_C = np.empty_like(t)
# record initial conditions
U_1[0] = z_0[0]
V_1[0] = z_0[1]
U_2[0] = z_0[2]
V_2[0] = z_0[3]
U_C[0] = z_0[4]
# solve ODE
for i in range(1,n):
    # span for next time step
    tspan = [t[i-1],t[i]]
    # solve for next step
    z = odeint(model,z_0,tspan,args=(varepsilon[i],))
    # store solution for plotting
    U_1[i] = z[1][0]
    V_1[i] = z[1][1]
    U_2[i] = z[1][2]
    V_2[i] = z[1][3]
    U_C[i] = z[1][4]
    # next initial condition
    z_0 = z[1]
B_1=(U_1+U_2)/2
B_2=(U_1-U_2)/2
# plot results
fin=int(finite_time/dist*n)
# plt.plot(t[:fin],B_1[:fin],'g',label='dark mode')
# plt.plot(t[:fin],varepsilon[:fin],'r:',label='vareps')
# plt.plot(t[:fin],B_2[:fin],'b-',label='light mode')
plt.ylabel('values')
plt.xlabel('time')
plt.legend(loc='best')
# plt.show()
