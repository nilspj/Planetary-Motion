import numpy as np
import common as co
import matplotlib.pyplot as plt 

#1.07603730251e-08#
h = 0.0001 #4.848137E-6 #4.848137E-6 #2.0846988E-5 
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

w = np.array([1/Ra,0])

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
    du2 = (1+alpha*u1**2)*GMs/(L**2) - u1
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

def magvect(v):
  magnitude = np.sqrt(v[0]*v[0] + v[1]*v[1])
  return magnitude

def angle(v1,v2):
  ang = np.arccos(np.dot(v1,v2)/(magvect(v1)*magvect(v2)))
  #ang = ang*180/np.pi
  return ang

def rungekutta(h,w,alpha):
    N = int(3*np.pi/h)
    #print(N)
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
        RKxList[i] = np.cos(i*h)/w[0]
        RKyList[i] = np.sin(i*h)/w[0]
        RKrList[i] = 1/w[0]
        RKvList[i] = 0 #np.sqrt(GMs*w[0])
        RKEkList[i] = 0.5*GMs*w[0]
        RKEpList[i] = -GMs*w[0]
        RKEmekList[i] = RKEkList[i] + RKEpList[i]
        w1 = RK4(p,h,w,alpha)
        w=w1
    Slice = RKrList[int(N/2):int(5*N/6):1]
    rpInd = int(N/2)+np.argmax(Slice)
    print(rpInd)
    precession = angle([RKxList[0],RKyList[0]],[RKxList[rpInd],RKyList[rpInd]])
    return RKxList,RKyList,precession,rpInd



alpha = 1.1*10**-8
print((100/T)*rungekutta(h,w,0.001)[2]*206264.81) #!!1.07603730251!! #1.07607391499
2
N = int(20)
aVec = np.linspace(alpha*10000,alpha*100000,N)
degVec = np.zeros(N)
for i in range(N):
    print(i)
    degVec[i] = (100/T)*rungekutta(h,w,aVec[i])[2]*206264.81
    
poly1 = np.polyfit(aVec,degVec,1)
print(poly1)
alpha = (43)/poly1[0]
print(alpha)
print((100/T)*rungekutta(h,w,alpha)[2]*206264.81)

#plt.ticklabel_format(useOffset=False)
plt.plot(aVec,degVec/1000, label='Precession per century',color='red')
plt.ylabel(r"$\theta\,[\mathrm{arcsec \times 10^{3}}]$",fontsize=14)
plt.xlabel(r"$\alpha\,[\mathrm{1}]$",fontsize=14)
plt.grid()
plt.legend()
plt.xlim(np.min(aVec),np.max(aVec))
'''
#MercuryA = analyticalEllipse(h)
MercuryN = rungekutta(h,w,0.01)
#orb1 = np.linspace(0,int(2*np.pi/h),int(2*np.pi/h)+1, dtype='int64')
#orb2 = np.linspace(3*int(2*np.pi/h),4*int(2*np.pi/h),int(2*np.pi/h)+1, dtype='int64')

#plt.plot(MercuryA[0],MercuryA[1], label='Mercury analytical orbit')
#plt.plot(np.take(MercuryN[0],orb1),np.take(MercuryN[1],orb1), label='Mercury numerical orbit1')
#plt.plot(np.take(MercuryN[0],orb2),np.take(MercuryN[1],orb2), label='Mercury numerical orbit2')
plt.plot(MercuryN[0],MercuryN[1])
#plt.plot(MercuryN[0][0],MercuryN[1][0], 'yo',markersize=10)
#plt.plot(MercuryN[0][MercuryN[3]],MercuryN[1][MercuryN[3]], 'yo',markersize=10)
plt.grid()
plt.plot(0,0, 'yo',markersize=10)
plt.legend(loc=1)
'''