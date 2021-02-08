import numpy as np 
GMs = 4*(np.pi**2)

def dist(x,y):
    return np.sqrt(x**2 + y**2)

def f(w):
    xdot = w[2]
    ydot = w[3]
    ax = -1*GMs*w[0]/(dist(w[0],w[1]))**3
    ay = -1*GMs*w[1]/(dist(w[0],w[1]))**3
    return np.array([xdot,ydot,ax,ay])
  
def RK4(f,h,w):
    s1 = f(w)
    s2 = f(w + h*s1/2)
    s3 = f(w + h*s2/2)
    s4 = f(w + h*s3)
    return (w + h*(s1 + 2*s2 + 2*s3 + s4)/6)

def EC(h,w):
    v1x = w[2] + f(w)[2]*h
    x1 = w[0] + v1x*h
    v1y = w[3] + f(w)[3]*h
    y1 = w[1] + v1y*h
    return np.array([x1,y1,v1x,v1y])

def eulercromer(h,T,w):
    N = int(T/h)+1
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
        ECEkList[i] = 0.5*dist(w[2],w[3])**2
        ECEpList[i] = -GMs/dist(w[0],w[1])
        ECEmekList[i] = ECEkList[i] + ECEpList[i]
        w1 = EC(h,w)
        w = w1
    return ECxList,ECyList,ECvxList,ECvyList,ECEkList,ECEpList,ECEmekList

def rungekutta(h,T,w):
    N = int(T/h)+1
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
        RKEkList[i] = 0.5*dist(w[2],w[3])**2
        RKEpList[i] = -GMs/dist(w[0],w[1])
        RKEmekList[i] = RKEkList[i] + RKEpList[i]
        w1 = RK4(f,h,w)
        w=w1
    return RKxList,RKyList,RKvxList,RKvyList,RKEkList,RKEpList,RKEmekList

def simplePlanetRK4(r,T,h): #### planets with circular orbits 
    y0 = 0
    x0 = r
    V0 = 2*np.pi*r/T#AU/Yr
    v0x =0  #AU/Yr
    v0y = V0 #AU/Yr
    w = np.array([x0,y0,v0x,v0y])
    return rungekutta(h,T,w)

def simplePlanetEC(r,T,h): #### planets with circular orbits 
    y0 = 0
    x0 = r
    V0 = 2*np.pi*r/T#AU/Yr
    v0x =0  #AU/Yr
    v0y = V0 #AU/Yr
    w = np.array([x0,y0,v0x,v0y])
    return eulercromer(h,T,w)

############ Three body system ####
Ms = 333480
#Mj = 0.1074*10000
#Me = 1
MJ = 317.89
Me = MJ
Mj = MJ
GMj = 4*(np.pi**2)*(Mj/Ms)
GMe = 4*(np.pi**2)*(Me/Ms)

def threeBodyF(w):
    Exdot = w[2]
    Eydot = w[3]
    Eax = -1*GMs*w[0]/(dist(w[0],w[1]))**3 - GMj*(w[0]-w[4])/(dist(w[4]-w[0],w[5]-w[1]))**3
    Eay = -1*GMs*w[1]/(dist(w[0],w[1]))**3 - GMj*(w[1]-w[5])/(dist(w[4]-w[0],w[5]-w[1]))**3
    Mxdot = w[6]
    Mydot = w[7]
    Max = -1*GMs*w[4]/(dist(w[4],w[5]))**3 - GMe*(w[4]-w[0])/(dist(w[4]-w[0],w[5]-w[1]))**3
    May = -1*GMs*w[5]/(dist(w[4],w[5]))**3 - GMe*(w[5]-w[1])/(dist(w[4]-w[0],w[5]-w[1]))**3
    return np.array([Exdot,Eydot,Eax,Eay,Mxdot,Mydot,Max,May])
              
def threeBodyRungekutta(h,T,w):
    N = int(T/h)+1
    ERKxList = np.zeros(N)
    ERKyList = np.zeros(N)
    MRKxList = np.zeros(N)
    MRKyList = np.zeros(N)
    ERKvxList = np.zeros(N)
    ERKvyList = np.zeros(N)
    ERKEkList = np.zeros(N)
    ERKEpList = np.zeros(N)
    ERKEmekList = np.zeros(N)
    MRKvxList = np.zeros(N)
    MRKvyList = np.zeros(N)
    MRKEkList = np.zeros(N)
    MRKEpList = np.zeros(N)
    MRKEmekList = np.zeros(N)
    for i in range(N):
        #t = t0 + h
        #t0 = t
        #r = dist(x0,y0)
        ERKxList[i]=w[0]
        ERKyList[i]=w[1]
        MRKxList[i]=w[4]
        MRKyList[i]=w[5]
        ERKvxList[i] = w[2]
        ERKvyList[i] = w[3]
        ERKEkList[i] = 0.5*dist(w[2],w[3])**2
        ERKEpList[i] = -GMs/dist(w[0],w[1])
        ERKEmekList[i] = ERKEkList[i] + ERKEpList[i]
        MRKvxList[i] = w[6]
        MRKvyList[i] = w[7]
        MRKEkList[i] = 0.5*dist(w[6],w[7])**2
        MRKEpList[i] = -GMs/dist(w[4],w[5])
        MRKEmekList[i] = MRKEkList[i] + MRKEpList[i]
        w1 = RK4(threeBodyF,h,w)
        w=w1
    return ERKxList,ERKyList,MRKxList,MRKyList

def threeBodySystemRK4(r1,T1,r2,T2,h,n): #### planets with circular orbits 
    Ey0 = 0
    Ex0 = r1
    EV0 = 2*np.pi*r1/T1#AU/Yr
    Ev0x =0  #AU/Yr
    Ev0y = EV0 #AU/Yr
    My0 = 0
    Mx0 = r2 
    MV0 = 2*np.pi*r2/T2#AU/Yr
    Mv0x =0  #AU/Yr
    Mv0y = MV0 #AU/Yr
    w = np.array([Ex0,Ey0,Ev0x,Ev0y,Mx0,My0,Mv0x,Mv0y])
    return threeBodyRungekutta(h,n*max(T1,T2),w)


