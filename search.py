import approx as apx
import numpy as np, matplotlib.pyplot as plt
import os
xi=1#np.sqrt(2)
def zfunc_nonlean_B(z, t, delta, alpha, gamma, gamma_C, s, g):
    b, B_20 = z
    # alpha = np.abs(xi)**2*( 3*chi_3/2-chi_2 ** 2/6)
    return [(-1j * delta- gamma ) * b + 1j * alpha*(np.abs(b) ** 2+2 * np.abs(B_20) ** 2) * b
            + 1j * alpha * B_20 ** 2 * np.conj(b),
            (-1j * delta- gamma)* B_20 + 1j * alpha * ((np.abs(B_20) ** 2 + 2 * np.abs(b) ** 2) * B_20
             + b ** 2 * np.conj(B_20)) + 1j / (1j + gamma_C) * (-s * B_20 + g / 2 / xi)]
B_bright=[]
B_dark=[]
delta=0.6#0.2581988897471611
alpha=1
gamma=0.1#0.11547005383792516
gamma_C=1e+02
s=1e-01
#39500
dist=500
phi=0
coef=5
start=0*coef#dist-100
stop= dist*coef
t = np.arange(0, dist,1/coef)
z0 = np.array([0.00001,0])#[4.441550009318063e-207, 0.8425760218073671])
g=36#24.921399999999977#15.495031461729917
for i in range(1400):
    g=g-0.01*i
    z, infodict = apx.odeintz(zfunc_nonlean_B, z0, t, args=(delta, alpha, gamma, gamma_C, s, g), full_output=True)
    B_1=z[:,0]*xi#*np.exp(1j*(1+delta)*t)
    B_2=z[:,1]*xi#*np.exp(1j*(1+delta)*t)
    print("B_1_stat = ",str(np.abs(B_1)[(dist-10)*coef]))
    print("B_2_stat = ",str(np.abs(B_2)[(dist-10)*coef]))
    print( r"g = " + str(g))
    B_bright.append((B_2)[(dist-10)*coef])
    B_dark.append((B_1)[(dist-10)*coef])
    if np.abs(B_dark[len(B_dark)-1])<2.0099314555365283e-192:
        z0 = np.array([2.0099314555365283e-192, B_bright[len(B_bright) - 1]])
    else:
        print(z0)
        z0 = np.array([np.real(B_dark)[len(B_dark)-1]+1j*np.imag(B_dark)[len(B_dark)-1], np.real(B_bright)[len(B_bright)-1]+1j*np.imag(B_bright)[len(B_bright)-1]])

    plt.clf()
    plt.plot(t[start:stop],(np.abs((B_2)))[start:stop], label='B_2_stat = '+ str(np.abs(B_2)[(dist-10)*coef]))
    plt.plot(t[start:stop], (np.abs((B_1)))[start:stop], label='B_1_stat = ' + str(np.abs(B_2)[(dist-10)*coef]))
    # plt.plot(t[start:stop],np.abs(B_2)[start:stop], label='B_2')
    # plt.plot(t[start:stop], np.abs(B_1)[start:stop], label='B_1')
    plt.text(500,0.2,r"g = " + str(g))
    print(max(np.abs(B_1)))
    plt.xlabel('t')
    plt.grid(True)
    plt.legend(loc='best')
    # lst=os.listdir(path="gisteresis2/down")
    a=[]
    plt.savefig("gisteresis2/"+str(i)+".png")
    with open('amplitudes.txt','a') as f:
        f.write("B_1_stat = " + str(np.abs(B_1)[(dist - 10) * coef])+  " B_2_stat = " + str(np.abs(B_2)[(dist - 10) * coef]) + r"g = " + str(g))
    # plt.show()



