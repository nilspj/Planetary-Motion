import numpy as np
import common as co
import matplotlib.pyplot as plt 

h = 0.001
N = int(2*np.pi/h)
Re = 6371000
GMe = 6.67*10**-11*6*10**24
Ms = 720 
Rp = 322000+Re #periphelion
Ra = Rp#35680000+Re # aphelion
e = (Ra-Rp)/(Ra+Rp)
a =Ra/(1+e)  #semi-major axis
b = a*np.sqrt(1-e**2)  # semi-minor axis
V0 = np.sqrt(GMe*(1+e)/(a*(1-e))) #velocity at periphelion(max)
V1= np.sqrt(GMe*(1-e)/(a*(1+e)))
L = Rp*V0 #angular momentum

w = np.array([1/Rp,0])

def radius(theta):
     r = (L**2)/(GMe*(1-e*np.cos(theta-np.pi)))
     return r

def g(w,L):
    u1 = w[0] # u = 1/r
    du1 = w[1]
    du2 = GMe/(L**2) - u1
    return np.array([du1,du2])


def RK4(f,h,w,L):
    s1 = f(w,L)
    s2 = f(w + h*s1/2,L)
    s3 = f(w + h*s2/2,L)
    s4 = f(w + h*s3,L)
    return (w + h*(s1 + 2*s2 + 2*s3 + s4)/6)

def analyticalEllipse(h):
    N = int(2*np.pi/h)
    xVec = np.zeros(N)
    yVec = np.zeros(N)
    for i in range(N):
        xVec[i] = radius(i*2*np.pi/N)*np.cos(np.pi + i*2*np.pi/N)
        yVec[i] = radius(i*2*np.pi/N)*np.sin(np.pi + i*2*np.pi/N)
    return xVec,yVec

def rungekutta(h,w,L):
    n = 1
    N = int(n*2*np.pi/h)
    RKxList = np.zeros(N)
    RKyList = np.zeros(N)
    RKrList = np.zeros(N)
    RKvList = np.zeros(N)
    RKEkList = np.zeros(N)
    RKEpList = np.zeros(N)
    RKEmekList = np.zeros(N)
    for i in range(N):
        
        RKxList[i]=np.cos(np.pi + n*i*2*np.pi/N)/w[0]
        RKyList[i]=np.sin(np.pi + n*i*2*np.pi/N)/w[0]
        RKrList[i] = 1/w[0]
        RKvList[i] = np.sqrt(GMe*(1+e)/(a*(1-e)))
        RKEkList[i] = 0.5*GMe*w[0]
        RKEpList[i] = -GMe*w[0]
        RKEmekList[i] = RKEkList[i] + RKEpList[i]
        w1 = RK4(g,h,w,L)
        w=w1
    
    return RKxList,RKyList

MercuryA = analyticalEllipse(h)
MercuryN = rungekutta(h,w,L)

#plt.plot(MercuryA[0],MercuryA[1], label='Mercury analytical orbit')
plt.plot(MercuryN[0]/1000000,MercuryN[1]/1000000,'-', label='Satellite numerical orbit', color='red')
plt.plot(0,0, 'bo',markersize=15, label='Earth')
plt.xlabel(r"$x \,[\mathrm{Mm}]$", fontsize=14)
plt.ylabel(r"$y \,[\mathrm{Mm}]$", fontsize=14)
plt.grid()
plt.legend(loc=1)