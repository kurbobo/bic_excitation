import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
finite_time=500#should be less then 50
delta_0=0.5
chi_21=1
chi_31=4/9
chi_22=chi_21
chi_32=chi_31
R=120
r_1=500
r_2=r_1
R_1=1
R_2=R_1
C_1=2E-7
C_2=C_1
C=1.1E-7
L_1=0.5E-3
L_2=L_1
zed=np.sqrt(L_1/C_1)
# gamma_U_1=np.sqrt(L_1/C_1)/r_1
varepsilon_0=0.3140916808149406
# k=1/np.sqrt(L_1*C_1)
# sigma=zed/R
# varepsilon_0=p/sigma/np.sqrt(2)
# print('zed=',zed)
# print('k=',k)
# alpha=3*chi_31/4+chi_21**2/6
# print('alpha=',alpha)
def model(z,t,varepsilon):
    noise =0# random.random()/10.0
    U_1 = z[0]+noise
    V_1 = z[1]
    U_2 = z[2]+noise
    V_2 = z[3]
    U_C = z[4]
    V_C = zed / R*(U_1 - U_2-U_C- varepsilon)
    dU_1dt = (1 + chi_21 * U_1 + chi_31 * U_1** 2)*(-zed / r_1 * U_1 + V_1 - V_C)
    dV_1dt = (-R_1 / zed * z[1] - z[0])
    dU_2dt = (1 + chi_22 * z[2] + chi_32 * z[2] ** 2)*(-zed * C_1 / r_2 / C_2 * z[2] + C_1 / C_2 * V_C - C_1 / C_2 * z[3])
    dV_2dt = (-R_2 * L_1 / zed / L_2 * z[3] + L_1 / L_2 * z[2])
    dU_Cdt =zed*C_1*(-1 / R / C * z[4] + 1 / C / R * (z[0] - z[2]) - 1 / C / R * varepsilon)
    dzdt = [dU_1dt, dV_1dt, dU_2dt, dV_2dt, dU_Cdt]
    return dzdt

# initial condition
z_0 = [0.01,0,0.0,0,0]

# number of time points
n = 500

# time points
dist=500
t = np.linspace(0,dist,n)

# step input
varepsilon = varepsilon_0*np.cos((1+np.cos(delta_0*t))*t)#*np.heaviside(20-t,0)

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
plt.plot(t[:fin],B_1[:fin],'g',label='dark mode')
# plt.plot(t[:fin],varepsilon[:fin],'r:',label='vareps')
plt.plot(t[:fin],B_2[:fin],'b-',label='light mode')
plt.ylabel('values')
plt.xlabel('time')
plt.legend(loc='best')
plt.show()
Gamma= 0.060000000000000005
alpha= 0.5
p= 0.18508029784877095
sigma= 0.13236581096849476
sigma2= 0.10027712952158695
