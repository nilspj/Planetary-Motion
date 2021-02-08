import numpy as np
import common as co
import matplotlib.pyplot as plt 

h = 0.001
N = int(2*np.pi/h) 
GMs = 4*(np.pi**2)
r = 0.39 #avarage distance
T = 0.2408 # Period
e = 0.206 #eccentricty 
Ra = 0.4667 #aphelion
Rp = 0.3075 # perihelion
a =Ra/(1+e)  #semi-major axis
b = a*np.sqrt(1-e**2)  # semi-minor axis
V0 = np.sqrt(GMs*(1+e)/(a*(1-e))) #velocity at periphelion(max)
L = Rp*V0 #angular momentum/planet mass 

w = np.array([1/Rp,0])

def radius(theta):
     r = (L**2)/(GMs*(1-e*np.cos(theta-np.pi)))
     return r

def g(w): ##### without precession
    u1 = w[0] # u = 1/r
    du1 = w[1]
    du2 = GMs/(L**2) - u1
    return np.array([du1,du2])

def p(w,alpha): ### with einstein's correction
    u1 = w[0] # u = 1/r
    du1 = w[1]
    du2 = GMs/(L**2)*(1+alpha*u1**2) - u1
    return np.array([du1,du2])

def RK4(f,h,w,alpha):
    s1 = f(w,alpha)
    s2 = f(w + h*s1/2,alpha)
    s3 = f(w + h*s2/2,alpha)
    s4 = f(w + h*s3,alpha)
    return (w + h*(s1 + 2*s2 + 2*s3 + s4)/6)

def analyticalEllipse(h):
    N = int(2*np.pi/h)
    xVec = np.zeros(N)
    yVec = np.zeros(N)
    for i in range(N):
        xVec[i] = radius(i*2*np.pi/N)*np.cos(np.pi + i*2*np.pi/N)
        yVec[i] = radius(i*2*np.pi/N)*np.sin(np.pi + i*2*np.pi/N)
    return xVec,yVec

def rungekutta(h,w,alpha):
    n=415 # antall runder i løpet av 100 år
    N = int(n*2*np.pi/h)
    RKxList = np.zeros(N)
    RKyList = np.zeros(N)
    RKrList = np.zeros(N)
    RKvList = np.zeros(N)
    RKEkList = np.zeros(N)
    RKEpList = np.zeros(N)
    RKEmekList = np.zeros(N)
    tList = np.zeros(N)
    for i in range(N):
        tList[i] = T/N*i
        RKxList[i] = np.cos(np.pi + i*n*2*np.pi/N)/w[0]
        RKyList[i] = np.sin(np.pi + i*n*2*np.pi/N)/w[0]
        RKrList[i] = 1/w[0]
        RKvList[i] = np.sqrt(GMs*w[0])
        RKEkList[i] = 0.5*GMs*w[0]
        RKEpList[i] = -GMs*w[0]
        RKEmekList[i] = RKEkList[i] + RKEpList[i]
        w1 = RK4(p,h,w,alpha)
        w=w1
        #print(w1[0])
    rpInd = np.argwhere(((Rp+0.0000001>RKrList)) & ((Rp-0.0000001)<RKrList))
    print(rpInd)
    #print(rpInd)
    precession = np.arccos(np.dot([RKxList[0],RKyList[0]],[RKxList[rpInd[-1]],RKyList[rpInd[-1]]])/(np.linalg.norm([RKxList[0],RKyList[0]])*np.linalg.norm([RKxList[rpInd[-1]],RKyList[rpInd[-1]]])))
    return RKxList,RKyList,precession,rpInd[-1]
'''
print(rungekutta(h,w,0)[2])
alpha = 1.1*10**(-8)
N = int(20)
aVec = np.linspace(alpha*1000,alpha*10000,N)
degVec = np.zeros(N)
for i in range(N):
    print(i)
    degVec[i] = rungekutta(h,w,aVec[i])[2]*206265
plt.ticklabel_format(useOffset=False)
plt.plot(aVec,degVec, label='Precession degree')
plt.grid()
'''

#MercuryA = analyticalEllipse(h)
MercuryN = rungekutta(h,w,0.0001)#1.07603730251*10**-8) # 1.07603730251*10**-8
orb1 = np.linspace(0,int(2*np.pi/h),int(2*np.pi/h)+1, dtype='int64')
orb2 = np.linspace(414*int(2*np.pi/h),415*int(2*np.pi/h),int(2*np.pi/h)+1, dtype='int64')
print(MercuryN[2])
#plt.plot(MercuryA[0],MercuryA[1], label='Mercury analytical orbit')
plt.plot(np.take(MercuryN[0],orb1),np.take(MercuryN[1],orb1), label='Mercury numerical orbit before precession')
plt.plot(np.take(MercuryN[0],orb2),np.take(MercuryN[1],orb2), label='Mercury numerical orbit after precession')
plt.plot(MercuryN[0][0],MercuryN[1][0], 'ko',markersize=5)
plt.xlabel(r"$x \,[\mathrm{AU}]$", fontsize=14)
plt.ylabel(r"$y \,[\mathrm{AU}]$", fontsize=14)
plt.plot(MercuryN[0][MercuryN[3]],MercuryN[1][MercuryN[3]], 'ko',markersize=5)
plt.grid()
plt.plot(0,0, 'o',markersize=15,label='Sun',color = 'gold')
plt.legend(loc=1)
