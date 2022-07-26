# -*- coding: utf-8 -*-
"""
Flow over a flat plate is solved here. 
Marching is done using Newton Raphson technique
"""
#%% Initialization
import math
import numpy as np
import matplotlib.pyplot as plt 
deta=0.0001
N=7/deta
n=int(N)
f=np.zeros(n)
g=np.zeros(n)
h=np.zeros(n)
eta=np.zeros(n) #non dimensional distance from wall
g[0]=0
f[0]=0
h0=0.1
h1=0.2
ind=np.zeros(n)
y=np.zeros(n)
U=1 # velocity m/s
nu=0.00001 # kinematic viscosity
#%% Newton Raphson marching
while abs(g[n-1]-1)>0.000001 :
    for i in range(1,n):
        h[0]=h0 
        f[i]=g[i-1]*deta+f[i-1]
        g[i]=h[i-1]*deta+g[i-1]
        h[i]=h[i-1]-0.5*f[i-1]*h[i-1]*deta
    h0=h[0]
    gp=g[n-1]
    for i in range(1,n):
        h[0]=h1   
        f[i]=g[i-1]*deta+f[i-1]
        g[i]=h[i-1]*deta+g[i-1]
        h[i]=h[i-1]-0.5*f[i-1]*h[i-1]*deta
    h1=h[0]  
    gn=g[n-1]
    h2=h1-((gn-1)*(h1-h0))/(gn-gp)
    h0=h1
    h1=h2
#%% Results
print(h1)
for i in range(1,n):
    eta[i]=eta[i-1]+deta
#plt.plot(g,eta)
#plt.xlabel('Non-dimensional velocity')
#plt.ylabel('Non-dimensional distance')
#%% Plot boundary layer
j=0
for i in range(1,n):
    if g[i]>=0.99:
        ind[j]=i
        j=j+1
indexForBL=int(ind[0])
BLconst=eta[indexForBL]
for x in range(1,n):
    # dist = 10 m. So 1/n = 0.001
    y[x]=BLconst*math.sqrt((nu*x/U))
#plt.plot(eta,y)
fig, axs= plt.subplots(2)
axs[0].plot(g, eta)
axs[0].set(xlabel='Non-dimensional velocity', ylabel='Non-dimensional distance',title='Dimless velocity from wall')
axs[1].plot(eta, y)
axs[1].set(xlabel='Distance from leading edge', ylabel='Boundary layer thickness',title='Boundary layer')