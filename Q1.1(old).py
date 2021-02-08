import numpy as np
import matplotlib.pyplot as plt 

N = 100
h = 1/N
#N = int(1/h)+1
y0 = 0
x0 = 1 
t0 = 0
V0 = 2*np.pi#AU/Yr
v0x =0  #AU/Yr
v0y = V0 #AU/Yr
r = 1 # AU
GMs = 4*(np.pi**2) #Au**3/yr**2
Msol = 1.989E30
w = np.array([x0,y0,v0x,v0y])
PI = np.pi
G = 6.67426
ME = (4*PI)/(G*333480)
GMS = 4*PI**2

def ri(xi, yi):
    if np.linalg.norm([xi,yi]) == 0:
        return 1
    else:
        return np.linalg.norm([xi,yi])

def dist(x,y):
    return np.sqrt(x**2 + y**2)

def f(w):
    xdot = w[2]
    ydot = w[3]
    ax = -1*GMs*w[0]/(dist(w[0],w[1]))**3
    ay = -1*GMs*w[1]/(dist(w[0],w[1]))**3
    return np.array([xdot,ydot,ax,ay])
  
def RK4(h,w):
    s1 = f(w)
    s2 = f(w + h*s1/2)
    s3 = f(w + h*s2/2)
    s4 = f(w + h*s3)
    return (w + h*(s1 + 2*s2 + 2*s3 + s4)/6)


def EC(h,w):
    w1 = f(w)
    v1x = w[2] + h*w1[2]
    v1y = w[3] + h*w1[3]
    x1 = w[0] + v1x*h
    y1 = w[1] + v1y*h
    return np.array([x1,y1,v1x,v1y])

########### EC ##############

w = np.array([x0,y0,v0x,v0y])

def eulercromer(h,w,N):
    N = int(N) + 1
    #N = int(1/h)+1
    ECxList = np.zeros(N)
    ECyList = np.zeros(N)
    ECvxList = np.zeros(N)
    ECvyList = np.zeros(N)
    ECEkList = np.zeros(N)
    ECEpList = np.zeros(N)
    ECEmekList = np.zeros(N)
    for i in range(N):
        #t = t0 + h
        #t0 = t
        #r = dist(x0,y0)
        ECxList[i]=w[0]
        ECyList[i]=w[1]
        ECvxList[i] = w[2]
        ECvyList[i] = w[3]
        ECEkList[i] = 0.5*(w[2]**2 + w[3]**2)
        ECEpList[i] = -GMs/np.sqrt((w[0]**2+w[1]**2))
        ECEmekList[i] = ECEpList[i] + ECEkList[i]
        w = EC(h,w)
    return ECxList,ECyList,ECvxList,ECvyList,ECEkList,ECEpList,ECEmekList

############ Reset############

y0 = 0
x0 = 1 
t0 = 0
V0 = 2*np.pi#AU/Yr
v0x =0  #AU/Yr
v0y = V0 #AU/Yr
r = 1 # AU

##### RK4 ####################

w = np.array([x0,y0,v0x,v0y])

def rungekutta(h,w,N):
    #N = int(1/h)+1
    N = int(N) + 1
    RKxList = np.zeros(N)
    RKyList = np.zeros(N)
    RKvxList = np.zeros(N)
    RKvyList = np.zeros(N)
    RKEkList = np.zeros(N)
    RKEpList = np.zeros(N)
    RKEmekList = np.zeros(N)
    for i in range(N):
        #t = t0 + h
        #t0 = t
        #r = dist(x0,y0)
        RKxList[i]=w[0]
        RKyList[i]=w[1]
        RKvxList[i] = w[2]
        RKvyList[i] = w[3]
        #RKEkList[i] = 0.5*dist(w[2],w[3])**2
        RKEkList[i] = 0.5*(w[2]**2 + w[3]**2)
        RKEpList[i] = -GMs/dist(w[0],w[1])
        RKEmekList[i] = RKEkList[i] + RKEpList[i]
        w = RK4(h,w)
    return RKxList,RKyList,RKvxList,RKvyList,RKEkList,RKEpList,RKEmekList
    
####### Plot ##############
y0 = 0
x0 = 1 
t0 = 0
V0 = 2*np.pi#AU/Yr
v0x =0  #AU/Yr
v0y = V0 #AU/Yr
r = 1 # AU

t_liste = np.linspace(0,1,num=int(1/h)+1)
plt.plot(rungekutta(h,w,N)[0],rungekutta(h,w,N)[1], label='Runge-kutta')
plt.plot(eulercromer(h,w,N)[0],eulercromer(h,w,N)[1], label='Euler-Cromer')
#plt.plot(rungekutta(0.1,w)[0],rungekutta(0.1,w)[1], label='RK4')
#plt.plot(eulercromer(0.1,w)[0],eulercromer(0.1,w)[1], label='EC')
#plt.plot(t_liste,eulercromer(h,w)[4], label='EC')
plt.plot(0,0,'o',markersize=15, label='Sun',color='gold')
plt.grid()
plt.legend(loc=1)
plt.xlabel(r"$x \,[\mathrm{AU}]$", fontsize=14)
plt.ylabel(r"$y \,[\mathrm{AU}]$", fontsize=12)
plt.figure()
plt.show()
#plt.xlim(-1.5,1.5)
#plt.ylim(-1.5,1.5)

####### Plot mekanisk energi med hensyn på h ##########
N = np.linspace(10,200,200)
N1 = np.linspace(10,200,10)
tau = 1/N #np.linspace(1/8,1/30,100)
tau1 = 1/N1
#N = 1/tau
RKmekList = np.zeros(len(tau))
ECmekList = np.zeros(len(tau))
RKmekList1 = np.zeros(len(tau1))
ECmekList1 = np.zeros(len(tau1))
for i in range(len(tau)):
    w = np.array([x0,y0,v0x,v0y])
    RKmek = rungekutta(tau[i],w,N[i])[6]
    ECmek = eulercromer(tau[i],w,N[i])[6]
    RKmekList[i] = RKmek[0]-RKmek[-1]
    ECmekList[i] = np.abs(ECmek[0]-ECmek[-1])
    if i % 20 == 0:
        j = int(i/20)
        RKmekList1[j] = RKmek[0]-RKmek[-1]
        ECmekList1[j] = np.abs(ECmek[0]-ECmek[-1])

print(ECmekList)
print(ECmekList1)

'''
print(np.argmax(RKmekList))
print(rungekutta(tau[np.argmax(RKmekList)],w)[4])
print(eulercromer(tau[np.argmax(RKmekList)],w)[4])
print(tau[np.argmax(RKmekList)])

xlist, ylist, vxlist, vylist, t_list, r_list = iterer_euler(x0, y0, delta_t, tmin, tmax, vx0, vy0)
velocity, kin_energy, pot_energy, tot_energy = energy_velocity(vxlist, vylist, r_list)
    #plt.plot(t_list, kin_energy, label = "Kinetisk energi")
    #plt.plot(t_list, pot_energy, label = "Potensiell energi")
plt.plot(t_list, tot_energy, label = "Total energi")
plt.title("Energi sfa. tid, "+ r'$\tau$ = ' + str(delta_t) + " startposisjon = (" + str(x0) + "," + str(y0) + ")")
plt.legend(loc = 'best')
plt.grid()
plt.show()
'''

#plt.plot(tau,RKmekList,label='Runge-kutta')
#plt.plot(tau,ECmekList,label='Euler-Cromer')
plt.semilogy(N,RKmekList,label='Runge-kutta')
#plt.semilogy(N,ECmekList,label='Euler-Cromer')
#plt.semilogy(N1,RKmekList1,label='Runge-kutta1')
plt.semilogy(N1,ECmekList1,label='Euler-Cromer')
#plt.plot(tau1,RKmekList1,label='RKmek1')
#plt.plot(tau1,ECmekList1,label='Euler-Cromer')
plt.legend(loc=1)
#plt.xlim(np.max(tau),np.min(tau))
plt.xlim(np.min(N),np.max(N))
#plt.xlabel(r"$\tau \,[\mathrm{Yr}]$", fontsize=14)
plt.xlabel(r"$N \,[\mathrm{1}]$", fontsize=14)
plt.ylabel(r"$\Delta E_{tot} \,[\mathrm{AU^2/Yr^2}]$", fontsize=12)
plt.figure()
plt.show()


    