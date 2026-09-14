
############################# Import Libraries ###############################
import numpy as np
import matplotlib.pyplot as plt
 
## Time course
dt = 0.1 #d
T = 100 #d
## Parameters - suppose Doxorubicin
ksc = 0.051*dt #/day stem proliferation rate
kcc = 0.105*dt #/day cancer proliferation rate
kwbc = 0.04*dt #/day WBC differentiation rate
kdsc = 0.01*dt #/day sc death rate
kdcc = 0.1*dt#/day cc death rate
kdwbc = 5*dt #/day wbc death rate
Awbc = 10**5 #WBC/SC
Dose = (110/1000)/(5*1000) #g/mL
kdecay = (-np.log(0.5)/0.6)*dt #/day
kconv = 6e16*dt#g/mL/d
kmax = 2e-16
kM = 10**-5#g/mL
## Initial Conditions
Nsc = [None] * int(T//dt)
Ncc = [None] * int(T//dt)
Nwbc = [None] * int(T//dt)
Cche = [None] * int(T//dt)
Nsc[1] = 10**6 #scs
Ncc[1] = 10**9 #ccs
Nwbc[1] = 8*10**8 #wbcs
Cche[1] = 0 #g/mL

## Setting Injection
InjFreq = np.array([10, 30, 50,70,90])/dt #d


## Run Simulation
for t in range(1,int(T//dt-1)):
    dNsc = ksc*Nsc[t]*(1 - kconv*(kmax*Cche[t])/(kM+Cche[t])) - Nsc[t]*kwbc - Nsc[t]*kdsc
    dNcc = kcc*Ncc[t]*(1 -kconv*(kmax*Cche[t])/(kM+Cche[t])) -Ncc[t]*kdcc
    dNwbc = Awbc*Nsc[t]*kwbc - Nwbc[t]*kdwbc
    dCche = -(Nsc[t]+Ncc[t])*(kmax*Cche[t])/(kM+Cche[t]) -Cche[t]*kdecay
    
    Nsc[t+1] = Nsc[t] + dNsc
    Ncc[t+1] = Ncc[t] + dNcc
    Nwbc[t+1] = Nwbc[t] + dNwbc
    Cche[t+1] = Cche[t] + dCche
    
    if sum(np.equal(t,InjFreq)):
        Cche[t+1] = Dose
    
## Plot Figure
plt.figure()
plt.subplot(411)
plt.plot(range(0,int(T//dt)), Cche)
plt.ylabel('Dox (g/mL)')
plt.subplot(412)
plt.plot(range(0,int(T//dt)), Nsc)
plt.ylabel('SCs')
plt.subplot(413)
plt.plot(range(0,int(T//dt)), Ncc)
plt.ylabel('Cancer')
plt.subplot(414)
plt.plot(range(0,int(T//dt)), Nwbc)
plt.ylabel('WBCs')
plt.xlabel('Time (d)')
plt.show()

