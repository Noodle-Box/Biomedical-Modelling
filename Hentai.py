import numpy as np
import matplotlib.pyplot as plt

""" n    
program prepared for BIOE3001 code demonstration 17/09/2024 by Seijun Stokes. 

----------------------make sure it can be run on Spyder!---------------------

Simulation 1- A(LV_DT): amount of Levosimendan in Digestive tract (mg)

Simulation 2- C(LV_SC), C(LV_OR): concentration of Levosimendan and active metabolite 
                                    OR-1896 in systemic circulation (mg/L)

Simulation 3- r(TM): radius of cerebral artery tunica media (mm)

Simulation 4- MAP: Mean Arterial Pressure (mmHg)

Simulation 5- v: cerebral blood flow velocity (cm/s)
"""

# define constants 

# dosing information
D_list = [1,4,8] # dosage of oral Levosimendan (mg) from (Nieminen et al., 2008) (Nieminen et al., 2013)
tau = 24 # dosing interval (h)

# Simulation 1 - digestive tract
V_SC = 3.75 # volume of systemic circulation (L)
KLV_DISS = 0.5 # Levosimendan dissolution constant (1/h)
KLV_ABS = 0.4 # Levosimendan absorbance constant (1/h)
t_s = 1 # start time of pill dissolution (h)
t_e = 3 # end time of pill dissolution (h)

# Simulation 2 - systemic circulation
t_LV_half = 1 # half-life of active Levosimendan (h)
KMET = 0.1 # Levosimendan metabolisation constant at gut (1/h)
t_OR_half = 80 # half-life of active metabolite OR-1896 (h)

# Simulation 3 - Middle Cerebral Artery tunica media
KLV_DIFF_TM = 0.05 # Levosimendan diffusion constant for MCA tunica media (1/h)
KLV_UPT_TM = 0.2 # Levosimendan uptake constant for MCA tunica media (1/h)
KOR_DIFF_TM = 0.05 # OR-1896 diffusion constant for MCA tunica media (1/h)
KOR_UPT_TM = 0.2 # OR-1896 uptake constant for MCA tunica media (1/h)
KDR_TM = 0.3 # Active drug agent multiplier constant in MCA tunica media (L*mm)/(mg*s)
gamma = 0.008 # radius ceiling function damping constant in tunica media (1/h)
r0 = 1.3 # tunica media radius initial value from (Kaspera et al., 2020; Mouches & Forkert, 2019) (mm)

# Simulation 4 - Myocardium
KLV_DIFF_CM = 0.05 # Levosimendan diffusion constant for cardiomyocytes (1/h)
KLV_UPT_CM = 0.01 # Levosimendan uptake constant for cardiomyocytes (1/h)
KOR_DIFF_CM = 0.05 # OR-1896 diffusion constant for cardiomyocytes (1/h)
KOR_UPT_CM = 0.01 # OR-1896 uptake constant for cardiomyocytes (1/h)
MAP_BL = 89 # Baseline (low bound) MAP.value taken from averaged S0157 (no stroke with hypertension)
MAP_0 = 114.0 # initial condition (high bound) MAP. value taken from averaged S0175 (stroke with hypertension)
KMAP = 150 # Sigmoidal curve multiplier constant (L/mg) 
C_MID = 0.04 # Midpoint (steepest drop) of the sigmoidal curve (mg/L)

# Simulation 5 - Middle Cerebral Artery Blood Flow velocity
ICP = 2666.44 # (20mmHg) assume constant intercranial pressure from (Jeon et al., 2014; Yuh & Dillon, 2010) (Pa)
L_CA = 0.025 # (25mm) assumed typical length of MCA. (Gunnal et al., 2019) (m)
nyu = 0.005 # (5 cP) blood viscosity. assume constant shear rate. from (Nader et al., 2019) (Pa*s)
KCORR = 0.025 # correction constant for calculating CBFV
v_0 = 44.0 # initial value of CBFv. value taken from averaged S0175 (stroke with hypertension)

# define delta function (approximation)
def delta(t:int, n:int, tau:int)->int:
    """
    approximation of dirac delta function in the context of drug delivery.
    returns 1 when (t-n*tau-t_s) = 0, returns 0 otherwise.
    """
    if t - n*tau - t_s == 0:
        return 1
    else:
        return 0
    
class DosageSimulator():
    """
    calculate values of various variables and plot them.
        
    reason for class object - ability to store values and compare later between instances
    """

    def __init__(self, num_days):
        """
        initialise parameters and variables
        """
        # number of hours to run the simulation
        self.num_hours = num_days*24

        # digestive tract
        self.A_LV_DT =[0.0] # amount of Levosimendan in digestive tract (mg)


        # systemic circulation
        self.C_LV_SC = [0.0] # concentration of Levosimendan in systemic circulation (mg/L)
        self.C_OR_SC = [0.0] # concentration of OR-1896 in systemic circulation (mg/L)

        # cerebral artery tunica media
        self.C_LV_TM = [0.0] # concentration of Levosimendan in tunica media (mg/L)
        self.C_OR_TM = [0.0] # concentration of OR-1896 in tunica media (mg/L)
        self.C_DR_TM = [0.0] # concentration of active agent (LV and OR) in tunica media (mg/L)
        self.r_TM = [r0] # radius of cerebral artery (mm)

        # myocardium
        self.C_LV_CM = [0.0] # concentration of Levosimendan in myocardium (mg/L)
        self.C_OR_CM = [0.0] # concentration of OR-1896 in myocardium (mg/L)
        self.C_DR_CM = [0.0] # concentration of active agent (LV and OR) in myocardium (mg/L)
        self.MAP = [MAP_0] # Mean Arterial Pressure (mmHg). initial value taken from S0175 (stroke with hypertension)

        # Middle Cerebral Artery blood flow
        self.v = [v_0] # bloodflow velocity in middle cerebral artery (cm/s). 
   
    def run_simulations(self, dosage_index:int):
        """
        simulate the change in each parameter according to time
        """
        D = D_list[dosage_index] # set the dosage amount for this unique simulation


        t = 0 # start the timer
        n = -1 # dosage number: dosing starts at 0

        while (t <= self.num_hours): 

            # Simulation 1 - Digestive tract
             
            if t % tau == 0:
                n += 1 # increment tau to the next dose

            """
            dA = Pill dosage + Pill dissolution - absorbtion through gut lumen
            """
            dA_LV_DT = KLV_DISS*(D*delta(t,n,tau)+np.heaviside((t-n*tau-t_s),0.5)*np.heaviside((t_e-t-n*tau),0.5))- KLV_ABS * (self.A_LV_DT[t] - V_SC * self.C_LV_SC[t])
            
            
            # Simulation 2 - Systemic circulation

            """
            dC_LV = absorbtion through intestinal lumen - diffusion to cardiomyocytes-diffusion to tunica media - excretion through liver - metabolisation at intestine

            dC_OR = metabolisation at intestine - diffusion to cardiomyocytes - diffusion to tunica media - excretion through urine
            """
            dC_LV_SC = KLV_ABS*(self.A_LV_DT[t]-V_SC*self.C_LV_SC[t])/V_SC - KLV_DIFF_CM*(self.C_LV_SC[t]-self.C_LV_CM[t]) - KLV_DIFF_TM*(self.C_LV_SC[t]-self.C_LV_TM[t]) - np.log(2)*self.C_LV_SC[t]/t_LV_half - KMET*self.C_LV_SC[t]

            dC_OR_SC = KMET*self.C_LV_SC[t] - KOR_DIFF_CM*(self.C_OR_SC[t]-self.C_OR_CM[t])- KOR_DIFF_TM*(self.C_OR_SC[t]-self.C_OR_TM[t]) - np.log(2)*self.C_OR_SC[t]/t_OR_half

            # Simulation 3 - MCA tunica media
            """
            dC =diffusion into TM - uptake by smooth muscle cell
            dr =  proportional to C_DR with ceiling
            """

            dC_LV_TM = KLV_DIFF_TM*(self.C_LV_SC[t]-self.C_LV_TM[t]) - KLV_UPT_TM*self.C_LV_TM[t]

            dC_OR_TM = KOR_DIFF_TM*(self.C_OR_SC[t]-self.C_OR_TM[t]) - KOR_UPT_TM*self.C_OR_TM[t]

            self.C_DR_TM.append(self.C_LV_TM[t] + self.C_OR_TM[t]) # combine the two active agents
            

            dr_TM = KDR_TM*self.C_DR_TM[t] - gamma*abs(self.r_TM[t] - r0)

            # Simulation 4 - Myocardium 
            """
            dC = diffusion into myocardium - uptake by cardiomyocytes
            MAP = Sigmoidal curve that starts at high value and gradually decreases to baseline.
            """
            dC_LV_CM = KLV_DIFF_CM*(self.C_LV_SC[t]-self.C_LV_CM[t]) - KLV_UPT_CM*self.C_LV_CM[t]

            dC_OR_CM = KOR_DIFF_CM*(self.C_OR_SC[t]-self.C_OR_CM[t]) - KOR_UPT_CM*self.C_OR_CM[t]

            self.C_DR_CM.append(self.C_LV_CM[t] + self.C_OR_CM[t]) # combine the two active agents
            

            self.MAP.append(MAP_BL + (MAP_0 - MAP_BL)/(1+ np.exp(KMAP*(self.C_DR_CM[t]-C_MID))))

            # Simulation 5 - Calculate the MCABFv
            """
            Use Poiseuille's law and fluid mechanics laws to find (deltaP)/(R*area)
            """
            self.v.append(KCORR *100 * (self.r_TM[t]*0.001)**2 * (self.MAP[t]*133.322 -20*0.77*133.322 - ICP) / (8 * nyu * L_CA))
            

            # Append all differentials to the value lists

            # Simulation 1
            self.A_LV_DT.append(self.A_LV_DT[t] + dA_LV_DT)

            # Simulation 2
            self.C_LV_SC.append(self.C_LV_SC[t] + dC_LV_SC)
            self.C_OR_SC.append(self.C_OR_SC[t] + dC_OR_SC)

            # Simulation 3
            self.C_LV_TM.append(self.C_LV_TM[t] + dC_LV_TM)
            self.C_OR_TM.append(self.C_OR_TM[t] + dC_OR_TM)
            self.r_TM.append(self.r_TM[t] + dr_TM)
            
            # Simulation 4
            self.C_LV_CM.append(self.C_LV_CM[t] + dC_LV_CM)
            self.C_OR_CM.append(self.C_OR_CM[t] + dC_OR_CM)


            # Progress in time by 1 hour
            t += 1
        
    def reset_values(self):
        """
        reset the key lists 
        """
        # digestive tract
        self.A_LV_DT =[0.0] # amount of Levosimendan in digestive tract (mg)


        # systemic circulation
        self.C_LV_SC = [0.0] # concentration of Levosimendan in systemic circulation (mg/L)
        self.C_OR_SC = [0.0] # concentration of OR-1896 in systemic circulation (mg/L)

        # cerebral artery tunica media
        self.C_LV_TM = [0.0] # concentration of Levosimendan in tunica media (mg/L)
        self.C_OR_TM = [0.0] # concentration of OR-1896 in tunica media (mg/L)
        self.C_DR_TM = [0.0] # concentration of active agent (LV and OR) in tunica media (mg/L)
        self.r_TM = [1.3] # radius of cerebral artery (mm)

        # myocardium
        self.C_LV_CM = [0.0] # concentration of Levosimendan in myocardium (mg/L)
        self.C_OR_CM = [0.0] # concentration of OR-1896 in myocardium (mg/L)
        self.C_DR_CM = [0.0] # concentration of active agent (LV and OR) in myocardium (mg/L)
        self.MAP = [114.0] # Mean Arterial Pressure (mmHg). initial value taken from S0175 (stroke with hypertension)

        # Middle Cerebral Artery blood flow
        self.v = [44.0] # bloodflow velocity in middle cerebral artery (cm/s). initial value taken from S0175 (stroke with hypertension)
   
    def plot_simulation_1(self):
        """
        digestive tract Levosimendan amount
        """
        plt.figure()

        # plot each dosage 
        for i, Dosage in enumerate(D_list):
            self.run_simulations(i)
            plt.plot(range(len(self.A_LV_DT)), self.A_LV_DT, label=f'{Dosage} mg')
            self.reset_values()

        # show the dosage interval
        x_values = np.arange(0, self.num_hours+2, tau) 
        for x_val in x_values:
            plt.axvline(x=x_val, color='r', linestyle='--', linewidth=1)

        plt.xlabel('Time (hours)')
        plt.ylabel('Levosimendan (mg)')
        plt.title('Amount of Levosimendan in digestive tract')
        # Display the legend for the dosage values
        plt.legend(title='Dose of Levosimendan (mg)', loc='upper right')

        plt.show() 

    def plot_simulation_2(self):
        """
        Plot systemic circulation concentrations of Levosimendan and OR-1896 for each dosage.
        """
        plt.figure()
        
        # Define a dictionary for Dosage -> Color mapping
        dosage_colors = {1: '#1f77b4', 4: '#ff7f0e', 8: '#2ca02c' }

        # Plot Levosimendan and OR-1896 for each dosage
        for i, Dosage in enumerate(D_list):
            self.run_simulations(i)

            # Get the color based on the dosage
            color = dosage_colors.get(Dosage, 'black')

            # Plot Levosimendan with a solid line and OR-1896 with a dashed line
            plt.plot(range(len(self.C_LV_SC)), self.C_LV_SC, label=f"Levosimendan {Dosage} mg dosage", color=color, linestyle='-')
            plt.plot(range(len(self.C_OR_SC)), self.C_OR_SC, label=f"OR-1896 for {Dosage} mg dosage of Levosimendan", color=color, linestyle='--')
            
            self.reset_values()

        # Show dosage intervals with vertical lines
        x_values = np.arange(0, self.num_hours + 2, tau)
        for x_val in x_values:
            plt.axvline(x=x_val, color='r', linestyle='--', linewidth=0.5)

        # Add labels and title
        plt.xlabel('Time (hours)')
        plt.ylabel('Concentration (mg/L)')
        plt.title('Levosimendan and OR-1896 in Systemic Circulation')

        # Display the legend for the dosage values, grouping Levosimendan and OR-1896
        plt.legend(title='Concentration of Levosimendan and OR-1896 (mg/L)', loc='upper right')

        # Show the plot
        plt.show()
  
    def plot_simulation_3(self):
        """
        MCA tunica media radius
        """
        plt.figure()

        # plot each dosage 
        for i, Dosage in enumerate(D_list):
            self.run_simulations(i)
            plt.plot(range(len(self.r_TM)), self.r_TM, label=f'{Dosage} mg')
            self.reset_values()

        # show the dosage interval
        x_values = np.arange(0, self.num_hours+2, tau) 
        for x_val in x_values:
            plt.axvline(x=x_val, color='r', linestyle='--', linewidth=1)

        plt.xlabel('Time (hours)')
        plt.ylabel('Tunica media radius (mm)')
        plt.title('MCA tunica media radius')
        # Display the legend for the dosage values
        plt.legend(title='Dose of Levosimendan (mg)', loc='upper right')

        plt.show() 

    def plot_simulation_4(self):
        """
        MAP
        """
        plt.figure()

        # plot each dosage 
        for i, Dosage in enumerate(D_list):
            self.run_simulations(i)
            plt.plot(range(len(self.MAP)), self.MAP, label=f'{Dosage} mg')
            self.reset_values()

        # show the dosage interval
        x_values = np.arange(0, self.num_hours+2, tau) 
        for x_val in x_values:
            plt.axvline(x=x_val, color='r', linestyle='--', linewidth=1)

        plt.xlabel('Time (hours)')
        plt.ylabel('MAP (mmHg)')
        plt.title('Mean Arterial Pressure')
        # Display the legend for the dosage values
        plt.legend(title='Dose of Levosimendan (mg)', loc='upper right')

        plt.show()

        

    def plot_simulation_5(self):
        """
        CBFv
        """
        plt.figure()

        # plot each dosage 
        for i, Dosage in enumerate(D_list):
            self.run_simulations(i)
            plt.plot(range(len(self.v)), self.v, label=f'{Dosage} mg')
            self.reset_values()

        # show the dosage interval
        x_values = np.arange(0, self.num_hours+2, tau) 
        for x_val in x_values:
            plt.axvline(x=x_val, color='r', linestyle='--', linewidth=1)

        plt.xlabel('Time (hours)')
        plt.ylabel('CBFv (cm/s)')
        plt.title('Cerebral Blood Flow velocity')
        # Display the legend for the dosage values
        plt.legend(title='Dose of Levosimendan (mg)', loc='upper right')

        plt.show()



if __name__ == "__main__":
    sim = DosageSimulator(18)
    sim.plot_simulation_1()
    sim.plot_simulation_2()
    sim.plot_simulation_3()
    sim.plot_simulation_4()
    sim.plot_simulation_5()
