##############################################################################
#
# BIOE3001: Quantitative Methods in Biomedical Engineering
# Tevyn Vergara | 46421201
#
# Main Code for Levosimendan Model
#
############################# Import Libraries ################################
import numpy as np
import matplotlib.pyplot as plt

################################# Time Course #################################
## Time course
dt = 0.1                           # hrs
T = 300                            # hrs, time course

###################### Parameters: Patient & Dosing ###########################

BW = 80                            # kg, body weight average from [INSERT SOURCE]
dose = 100                         # ng/kg/min, dosing rate
T_inf = 24                         # hrs, infusion duration
C_peak = 100                       # ng/mL, target plasma conc. at end of infusion

# Infusion runs CONSTANT for the full 24 hrs, so levosimendan reaches a steady
# state rather than ramping. V is set so that plateau = C_peak.
#   at steady state  C_ss = (dose*BW*60/V) / k0_per_hr
V = dose*BW*60/((np.log(2)/0.9)*C_peak)   # mL, apparent distribution volume
R = (dose*BW*60/V)*dt                     # ng/mL/step, infusion as a conc. rate

################### Parameters: Drug (Pharmacokinetics) #######################

# Parameters: Levosimendan (Active)
k0 = (np.log(2)/0.9)*dt           # 1/Hrs, Half life decay form [INSERT SOURECE]
f0 = 0.03                          # Free Fraction

# Parameters: Gut Transit Compartment (time capture, no activity)
#   Levosimendan is excreted into the small intestine, reduced by gut bacteria,
#   then reabsorbed. The gut holds the drug, which is why OR-1855 does not
#   appear until well after the parent has cleared.
T_gut = 30                         # hrs, 
kg = (np.log(2)/17.4)*dt           # 1/hrs, gut release rate

# Parameters: OR-1855 Metabolite (Unactive)
k1 = (np.log(2)/61)*dt             # 1/Hrs, Half life decay, from [INSERT SOURCE]
f1 = 0.58                          # Free Fraction OR-1855
m1 = 0.00451                       # Formation Fraction (Levo --> OR-1855)
                                   # ** TO PARAMETERISE ** fitted to the metabolite
                                   # data. Optimise later against the PD

# Parameters: OR-1896 Metabolite (Active)
k2 = (np.log(2)/32)*dt             # 1/Hrs, Half life decay ** FITTED **
                                   # NB Puttonen 2008 measured 62 +/- 5 hrs in
                                   # healthy subjects - the fit wants about half
f2 = 0.61                          # Free Fraction OR-1896

# Parameters: Back-Conversion (OR-1855 <--> OR-1896)
#   NAT2 acetylates OR-1855 into OR-1896, and OR-1896 is deacetylated back again.
#   The reaction runs BOTH ways - confirmed by [INSERT SOURCE].
#   This exchange is what lets both metabolites stay elevated at the same time,
#   which a one-way chain cannot reproduce.
kac = (np.log(2)/2.0)*dt           # 1/Hrs, acetylation   OR-1855 --> OR-1896
kbk = (np.log(2)/2.8)*dt           # 1/Hrs, deacetylation OR-1896 --> OR-1855


################### Parameters: Body (Pharmacodynamics) #######################



######################### Initial Conditions ##################################

CLevo = [None] * int(T//dt)           # ng/mL, Levosimendan
CGut  = [None] * int(T//dt)           # ng/mL, Gut transit compartment
C1855 = [None] * int(T//dt)           # ng/mL, OR-1855
C1896 = [None] * int(T//dt)           # ng/mL, OR-1896

CLevo[0] = 0                          # ng/mL, drug naive at t = 0
CLevo[1] = 0                          # ng/mL
CGut[0]  = 0                          # ng/mL
CGut[1]  = 0                          # ng/mL
C1855[0] = 0                          # ng/mL
C1855[1] = 0                          # ng/mL
C1896[0] = 0                          # ng/mL
C1896[1] = 0                          # ng/mL

############################# Main Loop #######################################

for t in range(1,int(T//dt-1)):

    ### Levosimendan (ACTIVE) ###

    # Constant infusion for the first 24 hours, so infusion and elimination run
    # together and levosimendan plateaus. After 24 hours, decay only.
    if t*dt < T_inf:
        dCLevo = R - k0*CLevo[t]
    else:
        dCLevo = -k0*CLevo[t]
    CLevo[t+1] = CLevo[t] + dCLevo

    ### Gut Transit Compartment (time capture) ###

    # Fills with the m1 fraction of what Levosimendan loses, from t = 0 onwards
    # (the drug is being eliminated throughout the infusion, not just after it)
    # Holds it, then releases into OR-1855 from T_gut onwards
    if t*dt < T_gut:
        Rel = 0
    else:
        Rel = kg*CGut[t]

    dCGut = m1*k0*CLevo[t] - Rel
    CGut[t+1] = CGut[t] + dCGut

    ### OR-1855 (UNACTIVE Metabolite) ###

    # acetylated forward, deacetylated back, plus its own elimination
    dC1855 = Rel - k1*C1855[t] - kac*C1855[t] + kbk*C1896[t]
    C1855[t+1] = C1855[t] + dC1855

    ### OR-1896 (ACTIVE Metabolite) ###

    # formed by acetylation, lost back to OR-1855 and by its own elimination
    dC1896 = kac*C1855[t] - kbk*C1896[t] - k2*C1896[t]
    C1896[t+1] = C1896[t] + dC1896


############################## Plot Figure ####################################

time = np.arange(0,int(T//dt))*dt  # hrs, x-axis in real time not step number

## Clinical data, OR-1855 from [INSERT SOURCE]
t_data = [24,25,26,28,30,36,37,38,40,44,48,60,72,96,144,192,240]
C_data = [0,0,0,0,0.05,0.8,0.9,1.05,1.2,1.55,2.1,1.9,1.95,1.5,0.7,0.35,0.15]

## Clinical data, OR-1896 from [INSERT SOURCE]
D_data = [0,0,0,0,0.1,0.6,0.7,0.85,1.1,1.75,2.35,2.4,2.35,2.05,1.2,0.55,0.3]

plt.figure()
plt.subplot(411)
plt.plot(time, CLevo)
plt.ylabel('Levosimendan (ng/mL)')
plt.subplot(412)
plt.plot(time, CGut)
plt.ylabel('Gut (ng/mL)')
plt.subplot(413)
plt.plot(time, C1855)
plt.plot(t_data, C_data, 'o')
plt.ylabel('OR-1855 (ng/mL)')
plt.subplot(414)
plt.plot(time, C1896)
plt.plot(t_data, D_data, 'o')
plt.ylabel('OR-1896 (ng/mL)')
plt.xlabel('Time (hrs)')
plt.tight_layout()
plt.show()

## Zoomed view, so the plateau and the 0.9 h decay are both visible
plt.figure()
plt.plot(time, CLevo)
plt.ylabel('Levosimendan ZOOM IN (ng/mL)')
plt.xlabel('Time (hrs)')
plt.show()

## OR-1855 against the clinical data
plt.figure()
plt.plot(time, C1855)
plt.plot(t_data, C_data, 'o')
plt.ylabel('OR-1855 vs DATA (ng/mL)')
plt.xlabel('Time (hrs)')
plt.show()

## OR-1896 against the clinical data
plt.figure()
plt.plot(time, C1896)
plt.plot(t_data, D_data, 'o')
plt.ylabel('OR-1896 vs DATA (ng/mL)')
plt.xlabel('Time (hrs)')
plt.show()

########################## Loading Patient Data ###############################