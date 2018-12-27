import numpy as np
from sympy import S, Eq, solve
from sympy import symbols
# delta_0=0#0.5
# chi_21=0#0.21
# chi_31=0#0.004
# chi_22=chi_21
# chi_32=chi_31
# R=10
# r_1=500
# r_2=r_1
# R_1=0.001
# R_2=R_1
# C_1=2E-5
# C_2=C_1
# C=2E-5
# L_1=1#0.5E-4
# L_2=L_1
# zed=np.sqrt(L_1/C_1)
# # gamma_U_1=np.sqrt(L_1/C_1)/r_1
# varepsilon=0.3140916808149406
# # gamma_V_1=R_1/np.sqrt(L_1/C_1)
# gamma=(np.sqrt(L_1/C_1)/r_1+R_1/np.sqrt(L_1/C_1))/2
# alpha=3 *chi_31/4+chi_21**2/6
# # gamma_C=np.sqrt(L_1/C_1)*C_1/R/C
# p=np.sqrt(2)*np.sqrt(L_1/C_1)*varepsilon/R
# # S=zed/R
# sigma=np.sqrt(L_1/C_1)/R/2/(1+(np.sqrt(L_1/C_1)*C_1/R/C)**2)
# sigma2=(C_1/R**2/C)/2/(1/np.sqrt(L_1/C_1)**2+(C_1/R/C)**2)
# # delta=0.5
# # sigma=0.102
# sigma2=0.1
# gamma=0.06
# alpha=0.5
# P=0.185
# k,m,sig,sig2,Pi,e,z,r,R_,Gamma=symbols('k m sig sig2 Pi e z r R_ Gamma',positive=True)
Gamma=0.06
Pi=0.185
# sig=0.102
sig2=0.1
k=4
# z=0.1
# r=1
z,r,m,e,R_,sig=symbols('z r m e R_ sig',positive=True)
equations=[Eq(m*e,Pi/pow(2,0.5)),
           Eq(z/r+R_/z,2*Gamma),
           Eq(1/m+k*k*m,1/2/sig),
           Eq(k/m/m+1/k,1/2/sig2),Eq(z,0.1),Eq(r,1)]
sol=dict((solve(equations))[0])
print(sol)
