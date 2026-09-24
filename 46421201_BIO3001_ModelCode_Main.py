# NOTE FOR MARKER: Approved Extension for 22/09/2026, 17:00

# Disclaimer: Claude Sonnet 5 was used in this assignment for:
#
# 1. Layout tips
# 2. Syntax checking
# 3. Reorganisation of structure, parameters and functions
# 4. Trouble shooting and debugging
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

############################### Time Course ###################################
## Time course
dt = 0.1                           # hrs
T = 300                            # hrs, time course

###################### Parameters: Patient & Dosing ###########################

# Dosing regimen: 0.1 ug/kg/min (100 ng/kg/min) constant IV infusion for 24
# hrs, no loading bolus. From [1]. 

BW = 80                            # kg, sample body weight
dose = 100                         # ng/kg/min, dosing infusion rate           
T_inf = 24                         # hrs, infusion duration          
C_peak = 100                       # ng/mL, max plasma concentration of infusion
                 
# Infusion conversions for main loop
V = dose*BW*60                     # ng/kg/min, infusion rate whole body rate

# Pre-defined half life here to solve for steady-state constant R
V = V/((np.log(2)/0.9)*C_peak)

# Infusion Rate going in during 24-hour infusion regimen
R = (dose*BW*60/V)*dt              # ng/mL, added per step

################### Parameters: Drug (Pharmacokinetics) #######################

# Parameters: Levosimendan (Active)
k0 = (np.log(2)/0.9)*dt            # 1/Hrs, half-life decay. From [1]
f0 = 0.03                          # Free Fraction. From [4]

# Parameters: Gut Transit Compartment (Time Delay between parent and metabolites)
T_gut = 30                         # hrs, time when gut starts to release (FITTED FROM [1]) 
kg = (np.log(2)/17.4)*dt           # 1/hrs, gut release rate (FITTED)

# Parameters: OR-1855 Metabolite (Unactive)
k1 = (np.log(2)/61)*dt             # 1/Hrs, half-life decay. From [1]
f1 = 0.58                          # Free Fraction OR-1855. From [4]
m1 = 0.0045                        # Formation Fraction from Levo (FITTED)
                                
# Parameters: OR-1896 Metabolite (Active)
k2 = (np.log(2)/62)*dt             # 1/Hrs, half-life decay. From [1]
f2 = 0.61                          # Free Fraction OR-1896. From [4]

# Parameters: Back-Conversion (OR-1855 <--> OR-1896) through NAT2 enzyme
kac = (np.log(2)/1.8)*dt           # 1/Hrs, acetylation (FITTED)
kbk = (np.log(2)/2.35)*dt           # 1/Hrs, deacetylation (FITTED)
                                  
################### Parameters: Body (Pharmacodynamics) #######################

# PD is kept in SI (Pa, m^3, Pa.s/m^3, m/s, Hz). PK above stays in ng/mL/hrs

# Conversion macros
MMHG_TO_PA = 133.322               # Pa per mmHg
ML_TO_M3   = 1e-6                  # m^3 per mL

# Blood and artery parameters
RHO    = 1060.0                    # kg/m^3, blood density. From [5]
    
### Parameters: Heart Block (Ce_eff --> SV, HR, MAP)
A_ao   = 7.07e-4                   # m^2, aortic cross-section. From [5]

# Hill Kinetics for drug binding to affect stroke volume
SV0     = 70.0 * ML_TO_M3            # m^3, baseline SV. From [7]
dSVmax  = 0.167                    # fractional max SV rise. From [3]
EC50_SV = 0.05                     # ng/mL (FITTED)

# Hill Kinetics for drug binding to affect heart rate [6]
HR0     = 70.0/60.0                # Hz, baseline heart rate
Emax_HR = 0.171                    # fractional max HR rise. From [2][3]
EC50_HR = 0.78                     # ng/mL, direct Hill term (FITTED)

# Hill Kinetics for drug binding to dilating system "body" vasculature
Emax_R   = 0.317                   # fractional max R_sys drop. From [3]
EC50_sys = 0.05                    # ng/mL (FITTED)

# Derivation for baseline systemic resistance 
# (how hard blood is push through vessels)
# Resistance = pressure / flow, both baseline
MAP0 = 93.0 * MMHG_TO_PA           # Pa, Baseline MAP. From [8]
CO0    = SV0 * HR0                 # m^3/s, baseline cardiac output (DERIVED)
R_sys0 = MAP0/CO0                  # Pa.s/m^3, baseline systemic resistance (DERIVED)

### Parameters: Brain Block (MAP, Ce_eff, CO --> V_cereb)

ICP_Pa = 10.0 * MMHG_TO_PA         # Pa, intracranial pressure. From [9]
R_MCA  = 0.0015                    # m, MCA radius. From [5]
A_MCA  = np.pi * (R_MCA)**2        # m^2, MCA cross-section

# Derivation for baseline cerebral resistance, similar to R_sys
# Resistance = pressure / flow, both baseline)
v_MCA0 = 0.48                      # m/s, baseline MCA velocity. From [2]
CPP0 = MAP0 - ICP_Pa               # Pa, baseline cerebral perfusssion (DERIVED)
Q_MCA0 = v_MCA0 * A_MCA            # m^3/s, baseline cerebral flow (DERIVED)
R_cereb0 = CPP0/Q_MCA0             # Pa.s/m^3, baseline cerebral resistance (DERIVED)

# Hill kinetics for drug binding to cerebral vasculature resistance
Emax_cereb = 0.065                 # fractional max R_cereb drop. From [2]
EC50_cereb = 0.03                  # ng/mL (FITTED)

# Effect-Compartment Delay to represent time taken for drug to reach and act on tissue
# Without it, PDs crash at t = 30
keo_thalf = 6.0                    # hrs, effect-compartment half-life (FITTED)
k_eo = (np.log(2)/keo_thalf)*dt    # 1/hrs, pre-multiplied by dt

### Heart (Mass Transfer System) ###
#
# Takes lagged effective drug and runs three Hill Kinetic responses 
# System resistance, stroke volume and heart rate. 
# Used to output MAP and cardiac output --> Brain
#
def heart_block(Ce_val):
    R_sys_v = R_sys0 * (1 - Emax_R*Ce_val/(EC50_sys + Ce_val))
    SV_v    = SV0 * (1 + dSVmax*Ce_val/(EC50_SV + Ce_val))
    HR_v    = HR0 * (1 + Emax_HR*Ce_val/(EC50_HR + Ce_val))
    CO_v    = SV_v * HR_v
    MAP_v   = CO_v * R_sys_v
    return R_sys_v, SV_v, HR_v, CO_v, MAP_v

### Cerebral Vasculature (Mass Transfer System) ###
#
# Takes lagged effecive drug, MAP and cardiac output
# Combined for Bernoulli's equation to solve for blood velocity
# V_cereb arranged in quadratic response to produce final output
# Equations from Dale lecture [5]
def cerebral_block(Ce_val, MAP_v, CO_v):
    Cerebral_R = R_cereb0*(1 - Emax_cereb*Ce_val/(EC50_cereb + Ce_val))
    Aortic_V = CO_v/A_ao
    A_term = Cerebral_R*A_MCA
    Blood_V = (-A_term + np.sqrt(A_term**2 + 2*RHO*(MAP_v-ICP_Pa) + RHO**2*Aortic_V**2)) / RHO
    return Cerebral_R, Aortic_V, Blood_V


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

Ce      = [None] * int(T//dt)         # ng/mL, active (plasma) drug concentration
Ce_eff  = [None] * int(T//dt)         # ng/mL, lagged effect-site concentration
R_sys   = [None] * int(T//dt)         # Pa.s/m^3, systemic resistance
SV      = [None] * int(T//dt)         # m^3, stroke volume
HR      = [None] * int(T//dt)         # Hz, heart rate
CO      = [None] * int(T//dt)         # m^3/s, cardiac output
MAP     = [None] * int(T//dt)         # Pa, mean arterial pressure
R_cereb = [None] * int(T//dt)         # Pa.s/m^3, cerebral resistance
V_ao    = [None] * int(T//dt)         # m/s, aortic velocity
V_cereb = [None] * int(T//dt)         # m/s, cerebral blood velocity

MAP[0] = 0                            # Pa
MAP[1] = 0
HR[0] = 0
HR[1] = 0
V_cereb[0] = 0
V_cereb[1] = 0

############################# Main Loop #######################################

# t = 0 and t = 1, matching CLevo[0:2] -- Ce and Ce_eff are drug-naive, so
# both PD blocks sit at their nominal baseline
Ce[0] = 0                             # ng/mL
Ce[1] = 0                             # ng/mL
Ce_eff[0] = 0                         # ng/mL
Ce_eff[1] = 0                         # ng/mL
for idx in (0, 1):
    (R_sys[idx], SV[idx], HR[idx], CO[idx], MAP[idx]) = heart_block(Ce_eff[idx])
    (R_cereb[idx], V_ao[idx], V_cereb[idx]) = cerebral_block(Ce_eff[idx], MAP[idx], CO[idx])

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

    if t*dt < T_gut:
        Rel = 0
    else:
        Rel = kg*CGut[t]

    dCGut = m1*k0*CLevo[t] - Rel
    CGut[t+1] = CGut[t] + dCGut

    ### OR-1855 (UNACTIVE Metabolite) ###
    
    # Formed by gut metabolism, dacetylation, lost from acetylation (1896) and self decay
    dC1855 = Rel - k1*C1855[t] - kac*C1855[t] + kbk*C1896[t]
    C1855[t+1] = C1855[t] + dC1855

    ### OR-1896 (ACTIVE Metabolite) ###

    # Formed by acetylation, lost back to OR-1855 and self decay
    dC1896 = kac*C1855[t] - kbk*C1896[t] - k2*C1896[t]
    C1896[t+1] = C1896[t] + dC1896

    ### Pharmacodynamics: Ce --> Ce_eff 

    Ce[t+1] = f0*CLevo[t+1] + f2*C1896[t+1]

    # Effect-Compartment for shape
    dCe_eff = k_eo*(Ce[t] - Ce_eff[t])
    Ce_eff[t+1] = Ce_eff[t] + dCe_eff

    ### Heart (Mass Transfer System) ###
    (R_sys[t+1], SV[t+1], HR[t+1], CO[t+1], MAP[t+1]) = heart_block(Ce_eff[t+1])

    ### Cerebral Vasculature (Mass Transfer System) ###
    (R_cereb[t+1], V_ao[t+1], V_cereb[t+1]) = cerebral_block(Ce_eff[t+1], MAP[t+1], CO[t+1])


############################## Plot Figure ####################################

time = np.arange(0,int(T//dt))*dt  # hrs, x-axis in real time not step number

## Clinical data of parent drug and metabolites, from Table 1 data
t_data    = [24,25,26,28,30,36,37,38,40,44,48,60,72,96,144,192,240]
levo_data = [100, 50, 20, 5, 1.3, 0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
m1_data   = [0,0,0,0,0.05,0.8,0.9,1.05,1.2,1.55,2.1,1.9,1.95,1.5,0.7,0.35,0.15]
m2_data   = [0,0,0,0,0.1,0.6,0.7,0.85,1.1,1.75,2.35,2.4,2.35,2.05,1.2,0.55,0.3]

## Levosimendan against clinical data [1]
plt.figure()
plt.plot(time, CLevo)
plt.plot(t_data, levo_data, 'o')
plt.ylabel('Levosimendan ZOOM IN (ng/mL)')
plt.xlabel('Time (hrs)')
plt.title('Levosimendan (Parent) vs Clinical Data')
plt.show()

## OR-1855 against the clinical data [1]
plt.figure()
plt.plot(time, CGut)
plt.ylabel('Gut Concentration (ng/mL)')
plt.xlabel('Time (hrs)')
plt.title('Gut concentration of Levosimendan (Time Compartment')
plt.show()

## OR-1896 against the clinical data [1]
plt.figure()
plt.plot(time, C1896)
plt.plot(t_data, m2_data, 'o')
plt.ylabel('OR-1896 vs DATA (ng/mL)')
plt.xlabel('Time (hrs)')
plt.title('OR-1896 (Active Metabolite) vs Clinical Data')
plt.show()

## Pharmacodynamic outputs -- converted from SI to clinically readable units
## just for these plots (the model itself stays in SI throughout)
MAP_mmHg    = np.array(MAP)/MMHG_TO_PA
HR_bpm      = np.array(HR)*60
V_cereb_cms = np.array(V_cereb)*100

t_settled = 2.0                       # hrs
mask = time >= t_settled

plt.figure()
plt.plot(time[mask], MAP_mmHg[mask])
plt.ylabel('MAP (mmHg)')
plt.xlabel('Time (hrs)')
plt.title('Mean Arterial Pressure (MAP)')
plt.show()

plt.figure()
plt.plot(time, HR_bpm)
plt.ylabel('HR (bpm)')
plt.xlabel('Time (hrs)')
plt.title('Heart Rate (HR)')
plt.show()

plt.figure()
plt.plot(time[mask], V_cereb_cms[mask])
plt.ylabel('Blood Velocity (cm/s)')
plt.xlabel('Time (hrs)')
plt.title('Cerebral Blood Flow Velocity')
plt.show()

########################## Loading Patient Data ###############################

# TBA

############################## References #####################################
# Declaration of AI use
# [0] Anthropic, "Claude, Sonnet 5," Claude, 18 September 2026 [Online]. https://claude.ai/new 
#
# Provided Table sources from Mark project context
#
# Full references details are in reference manager 
#
# Table 1. 
# [1] {Puttonen, 2008}
#
# Table 2. 
# [2] {Kivikko, 2015}
#
# Table 3. 
# [3] {Álvarez, 2006}
#
# Others:
#
# [4] {Antila, 2004}
# [5] {Robinson, 2026}
# [6] {Kashou, 2022}
# [7] {Khalil, 2019}
# [8] {DeMers, 2023}
# [9] {Munakomi, 2026}
