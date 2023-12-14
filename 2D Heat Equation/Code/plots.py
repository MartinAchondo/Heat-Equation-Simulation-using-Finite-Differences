
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi'] = 200
plt.rcParams['savefig.dpi'] = 200


# Grafico de temperatura ---------

T = pd.read_csv("results/Results_T_eu.txt", header=None, delim_whitespace=True)
Lx = 2
Ly = 1
Ny,Nx = T.shape

x = np.linspace(0,Lx,Nx)
y = np.linspace(0,Ly,Ny)

X,Y = np.meshgrid(x,y)

fig0, ax0 = plt.subplots(1, 1, )
cf0 = ax0.contourf(X,Y,T,np.arange(0, 80, .1), cmap="hot")
cbar0 = plt.colorbar(cf0,)

cbar0.set_label('Temperatura [°C]')
fig = plt.gcf()
fig.set_size_inches(8, 4)

name = 'Temperature'
path_dir = os.path.join(os.getcwd(),name+'.png')
fig.savefig(path_dir)


# Grafico de flujo de calor ------

q = pd.read_csv("results/Results_q_eu.txt", header=None, delim_whitespace=True)
qx = pd.read_csv("results/Results_qx_eu.txt", header=None, delim_whitespace=True)
qy = pd.read_csv("results/Results_qy_eu.txt", header=None, delim_whitespace=True)

Lx = 2
Ly = 1
Ny,Nx = q.shape

x = np.linspace(0,Lx,Nx)
y = np.linspace(0,Ly,Ny)
X,Y = np.meshgrid(x,y)

plt.contourf(X,Y,q/1000,300, cmap="hot")
cbar = plt.colorbar()
plt.quiver(X,Y,qx*100,qy*100,color='k')
fig = plt.gcf()
cbar.set_label('Flujo calor [kW/m2]')
fig = plt.gcf()
fig.set_size_inches(8, 4)
plt.show()
name = 'heat'
plt.savefig(os.path.join(os.getcwd(),name+'.png'))


# Grafico de parte simetrica --------

T_b = np.zeros(Nx+Nx+Ny-2)
k = -1
Tx = T.to_numpy()

T_b[:Nx] = Tx[0,:]
T_b[Ny//2+Nx+1:Ny//2+Nx+Nx] = Tx[Ny//2,-1:0:-1]
kk = 0
T2 = np.zeros(Nx+Nx+Ny-2+Ny//2)
for k in range(len(T_b)-1,0,-1):
    T2[kk] = T_b[k]
    kk = kk+1

ll = np.linspace(0,5,253)
plt.plot(ll,T2[Ny//2-2:],linewidth=1,color='k', label='Temperatura')
plt.grid()
plt.legend(loc='upper left')
plt.xlabel('L [m]')
plt.ylabel('T [°C]')
plt.ylim([-1,80])
fig = plt.gcf()
fig.set_size_inches(6, 4)
plt.show()
name = 'tsym'
plt.savefig(os.path.join(os.getcwd(),name+'.png'))
