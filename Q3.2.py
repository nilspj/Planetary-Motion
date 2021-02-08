import numpy as np
import common as co
import matplotlib.pyplot as plt 

h = 0.001
N = int(2*np.pi/h)
Re = 6371000
GMe = 6.67*10**(-11)*5.972*10**(24)
Ms = 720
Rp = 322000+Re #periphelion
Ra = Rp # aphelion
e = (Ra-Rp)/(Ra+Rp)
a =Ra/(1+e)  #semi-major axis
b = a*np.sqrt(1-e**2)  # semi-minor axis
V0 = np.sqrt(GMe/Rp) #velocity at periphelion(max)
V1= np.sqrt(GMe*(1-e)/(a*(1+e)))
L = Rp*V0 #angular momentum


def radius(theta):
     r = (L**2)/(GMe*(1-e*np.cos(theta-np.pi)))
     return r

w = np.array([1/Rp,0]) # initial conditions for u and u'
def g(w,L):
    u1 = w[0] # u = 1/r
    du1 = w[1] # u'
    du2 = GMe/(L**2) - u1 # u''
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
    n = 3 ## number of orbits
    N = int(n*2*np.pi/h) 
    Rp = 322000+Re #periphelion
    Ra = Rp # aphelion
    RKxList = np.zeros(N)
    RKyList = np.zeros(N)
    RKrList = np.zeros(N)
    RKvList = np.zeros(N)
    RKEkList = np.zeros(N)
    RKEpList = np.zeros(N)
    RKEmekList = np.zeros(N)
    DeltaV = 0
    for i in range(N):
        e = (Ra-Rp)/(Ra+Rp)
        a = Ra/(1+e)  #semi-major axis
        RKxList[i]=np.cos(np.pi + n*i*2*np.pi/N)/w[0]
        RKyList[i]=np.sin(np.pi + n*i*2*np.pi/N)/w[0]
        RKrList[i] = 1/w[0]
        RKvList[i] = np.sqrt(GMe*(1-e)/(a*(1+e)))
        #print(RKvList[i])
        #print(w[0])
   
        RKEkList[i] = 0.5*GMe*w[0]
        RKEpList[i] = -GMe*w[0]
        RKEmekList[i] = RKEkList[i] + RKEpList[i]
        w1 = RK4(g,h,w,L)
        w=w1
        
        if i == int(N/3): ### first increase in velocity 
            Ra = 35680000+Re
            L= Rp*(RKvList[i]+np.sqrt(GMe/Rp)*(np.sqrt(2*Ra/(Ra+Rp))-1))
            DeltaV += L/Rp-RKvList[i]
        if i == int(0.5*N): ####### second increase in velocity
            #### if we don't multiply with a number lower than 1 the orbit diverges
            L = Ra*(RKvList[i]+np.sqrt(GMe/Ra)*(1-np.sqrt(2*Rp/(Ra+Rp)))) #0.676*
            DeltaV += L/Ra-RKvList[i]
            Rp = Ra

    tH = np.sqrt(np.pi**2*(35680000+Re+322000+Re)**3/(8*GMe))
    return RKxList,RKyList,tH,DeltaV

def rungekutta2(h,w,L):
    n = 6 ## number of orbits
    N = int(n*2*np.pi/h) 
    Rp = 322000+Re #periphelion
    Ra = Rp # aphelion
    RKxList = np.zeros(N)
    RKyList = np.zeros(N)
    RKrList = np.zeros(N)
    RKvList = np.zeros(N)
    RKEkList = np.zeros(N)
    RKEpList = np.zeros(N)
    RKEmekList = np.zeros(N)
    DeltaV = 0
    for i in range(N):
        e = (Ra-Rp)/(Ra+Rp)
        a = Ra/(1+e)  #semi-major axis
        RKxList[i]=np.cos(np.pi + n*i*2*np.pi/N)/w[0]
        RKyList[i]=np.sin(np.pi + n*i*2*np.pi/N)/w[0]
        RKrList[i] = 1/w[0]
        RKvList[i] = np.sqrt(GMe*(1-e)/(a*(1+e)))
        #print(RKvList[i])
        #print(w[0])
   
        RKEkList[i] = 0.5*GMe*w[0]
        RKEpList[i] = -GMe*w[0]
        RKEmekList[i] = RKEkList[i] + RKEpList[i]
        w1 = RK4(g,h,w,L)
        w=w1
        if i == int(N/6): ### first increase in velocity 
            Ra = 20000000+Re
            L= Rp*(RKvList[i]+np.sqrt(GMe/Rp)*(np.sqrt(2*Ra/(Ra+Rp))-1))
            DeltaV += L/Rp-RKvList[i]
        if i == int(N/4): ####### second increase in velocity
            #### if we don't multiply with a number lower than 1 the orbit diverges
            L = Ra*(RKvList[i]+np.sqrt(GMe/Ra)*(1-np.sqrt(2*Rp/(Ra+Rp)))) #0.676*
            DeltaV += L/Ra-RKvList[i]
            Rp = Ra
        
        if i == int(N/2): ### first increase in velocity 
            Ra = 35680000+Re
            L= Rp*(RKvList[i]+np.sqrt(GMe/Rp)*(np.sqrt(2*Ra/(Ra+Rp))-1))
            DeltaV += L/Rp-RKvList[i]

        if i == int(7*N/12): ####### second increase in velocity
            #### if we don't multiply with a number lower than 1 the orbit diverges
            L = Ra*(RKvList[i]+np.sqrt(GMe/Ra)*(1-np.sqrt(2*Rp/(Ra+Rp)))) #0.676*
            DeltaV += L/Ra-RKvList[i]
            Rp = Ra

    tH0 = np.sqrt(np.pi**2*(20000000+Re+322000+Re)**3/(8*GMe))
    tH1 = np.sqrt(np.pi**2*(35680000+Re+20000000+Re)**3/(8*GMe))
    tH = tH0 + tH1
    return RKxList,RKyList,tH,DeltaV


SatA = analyticalEllipse(h)
SatN = rungekutta(h,w,L)
SatN2 = rungekutta2(h,w,L)
print(SatN[2],'time')
print(SatN[3],'Delta v')
print(SatN2[2],'time new')
print(SatN2[3],'Delta v new')
#plt.plot(SatN2[0],SatN2[1],'g--', label='Satellite numerical orbit')
plt.plot(SatN[0],SatN[1],'r--', label='Satellite numerical orbit')
plt.plot(0,0, 'bo',markersize=15,label='Earth')
plt.xlabel(r"$x \,[\mathrm{m}]$", fontsize=14)
plt.ylabel(r"$y \,[\mathrm{m}]$", fontsize=14)
plt.grid()
plt.legend(loc=1)
plt.show()