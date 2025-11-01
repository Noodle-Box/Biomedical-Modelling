###################################################################################
# Project File for BIOE3001 Project
# Modelling the effect of Dofetilide on diagnosed Short QT patients
#
# Author: Tevyn Vergara
#
#
# ChatGPT used to restructure in class and other functionalities - based simulation.
# Reference:
# OpenAI. (2025). ChatGPT 4o. [Large language model].
# https://chat.openai.com/chat

"""
Simulation 1: C(Dof): concentration of Dofetilide in blood (mg/L)

Simulation 2: IKr(Block): Potassium Channel Current (IKr) reduction from to Dofetilide binding (uAmps)

Simulation 3: t(repolar): Total Action potential repolarization time (ms)

Simulation 4: t(QT): QT Interval duration (ms)
"""
###################################################################################
# Standard library imports

import numpy as np
import matplotlib.pyplot as plt
from typing import Optional 

###################################################################################
### Parameters for Simulation ###

# Three dosing regimens to compare in plots/printing
REGIMENS = [
    {"label": "250 mg x2 then 125 mg, q12h", "doses": (250.0, 125.0)},
    {"label": "125 mg x2 then 125 mg, q12h", "doses": (125.0, 125.0)},
    {"label": "100 mg x2 then 50 mg, q12h",  "doses": (100.0, 50.0)},
]

# -----------------------------
# Schedule: 250 mg ×2 then 125 mg, q12h
# -----------------------------
tau = 12                          # dosing interval (hrs) 
ts = 0                            # starting time for first (hrs)

# Heaviside dosing function
def delta(t: int, n: int, tau: int) -> int:
    return 1 if (t - n*tau - ts) == 0 else 0

# -----------------------------
# PK parameters (units in comments)
# -----------------------------
V   = 70.0                     # (L)      apparent central volume
kAbsorptionBS   = 0.6                      # (1/h)    first-order absorption from gut, from (...)
KDISS  = 0.8                      #          tablet absorption efficiency
t_half = 10.0                     # h        elimination half-life, from (...)
kDecay  = np.log(2) / t_half/2    # 1/h 

kAbsorption = (1.0 / 1.5)                    # (1/h)    absorption rate constant, from (...)

# Extra sink due to binding removing drug from plasma
KBinding = 0.0001                 # (1/h)    binding kinetic of Dofetillide to IKr channel

# Extra concentration-driven block term to increase regimen separation (dimensionless)
# Scales additional fractional block with concentration via a saturating function.
kBlock = 0.25                     #          fractional block at high concentration

# -----------------------------
# Binding (Eq. 1): Block(C) in [0,1]
# -----------------------------

""" Tuning pameters for parameterisation """
k1 = 5                           #           Emax (≤1.0), part of Kbinding calculation
k2 = (5e-3)
k3 = 3e2                         #           dimensionless gain on repolarization scaling
k4 = 1.5                         #           curvature exponent in Eq. 3 (tune)
k5 = 4e-1                      #           scaling from Δt_repolar to QT (ms per ms * factor in Eq. 4)

# Hill Mechanics Parameters
EC50 = 20                        # (mg/L)     — tune with your data
hill = 0.015                     # (n)       (Hill coefficient, from (...)

# Hill function for binding mechanism
def hill_mech(v, Emax = k1, EC50=EC50, n=hill):

    v = max(v, 1e-12)
    num = Emax * (v ** n)
    den = (EC50 ** n + v ** n)
    return float(num / den)

# -----------------------------
# Electrophysiology (Eq. 2–4)
# -----------------------------
IKR0_uA       = 200*10**-12     # (A)          baseline IKr magnitude, from (...)
t0_repolar_ms = 160.0           # (ms)         baseline repolarization time, from (...)
t0_QT_ms      = 310.0           # (ms)         baseline QT, from (...)
EPS = 1e-9

""" Main Model Simulation """
class DosageSimulator:
    """
    Discrete-hour simulator (mirrors your reference structure and the previous code),
    with states and outputs recorded each hour.
    """
    def __init__(self, dose_hours: int, regimen: tuple = (250.0, 125.0), label: Optional[str] = None):

        self.num_hours = dose_hours
        # Regimen: (first two doses, subsequent)
        self.set_regimen(regimen, label)
        self.initial_values()

    def set_regimen(self, regimen: tuple, label: Optional[str] = None):
        self.regimen_first, self.regimen_after = float(regimen[0]), float(regimen[1])
        self.regimen_label = label if label is not None else f"{int(self.regimen_first)} mg x2 then {int(self.regimen_after)} mg, q12h"

    def dose_at(self, n: int) -> float:
        return self.regimen_first if n < 2 else self.regimen_after

    def initial_values(self):

        # Blood drug concentration
        self.C_dof = [0.0]                  # (mg/L), in plasma

        # Potassium current 
        self.IKr = [IKR0_uA]                # (µA) , post-block current in Eq. 2

        # Cardiac AP repolarization time
        self.t_repolar = [t0_repolar_ms]    # (ms)

        # QT interval duration
        self.t_QT_ms = [t0_QT_ms]           # (ms)

    def run_simulations(self):
        t = 0
        n = -1  # dose window index

        while t <= self.num_hours:
            if t % tau == 0:
                n += 1 # increment tau to the next dose
            D = self.dose_at(n)  # mg scheduled for this window depending on regimen

            """ Simulation 1 - Dofetilide concentration in blood (mg/L) """

            # Input dose calculation with first-order absorption  
            Fin = (D / V) * kAbsorption * np.exp(-kAbsorption * (t % tau))
            C_new = self.C_dof[t]*(1 - kDecay - KBinding) + Fin
            self.C_dof.append(C_new)

        
            """ Simulation 2 - IKr after Dofetilide block (pAmps)"""

            # Concentration value and concentration-driven term to calculate IKr
            C_eff = self.C_dof[-1]
            r_bind = KBinding * C_eff
            block_base = (hill_mech(r_bind)/k1 - k2*kDecay*C_eff)
            # Additional saturating concentration term to help differentiate regimens
            conc_term = kBlock * (C_eff / (EC50 + C_eff)) 
            block = block_base + conc_term
            # Calculation of IKr after dofetilide bindng
            IKr_current = (IKR0_uA * (1.0 - block))
            self.IKr.append(IKr_current)

            """ Simulation 3 - Repolarization time of ventricular A.P (ms) """
            t_repolar = (t0_repolar_ms + (k3 * (block)**k4)) 
            self.t_repolar.append(t_repolar)

            """ Simulation 4 - QT Interval duration (ms) """
            t_QT = t0_QT_ms + k5 * (self.t_repolar[t] - t0_repolar_ms) 
            self.t_QT_ms.append(t_QT)

            # Increment time
            t += 1
        # -----------------------------
        # Plotters (your 4 simulations)
        # -----------------------------
    
    # Plot the concentration of Dofetilide in blood over time
    def plotsimulation_1(self):
        """ Simulation 1: C(Dof) in blood (mg/L) """
        plt.figure()
        for reg in REGIMENS:
            self.set_regimen(reg["doses"], reg["label"]) ; self.initial_values(); self.run_simulations()
            plt.plot(range(len(self.C_dof)), self.C_dof, label=self.regimen_label)
        for x in np.arange(0, self.num_hours + 2, tau):
            plt.axvline(x=x, color="r", linestyle="--", linewidth=0.7)
        plt.xlabel("Time (hours)"); plt.ylabel("C(Dof) (mg/L)")
        plt.title("Simulation 1: Dofetilide concentration in blood")
        plt.legend(); plt.show()
        self.initial_values()

    # Plot the IKr current after Dofetilide binding
    def plotsimulation_2(self):
        """Simulation 2: IKr after Dofetilide block (A)."""
        plt.figure()
        for reg in REGIMENS:
            self.set_regimen(reg["doses"], reg["label"]) ; self.initial_values(); self.run_simulations()
            series = self.IKr
            ylab = "IKr (A)"
            ttl = "Simulation 2: IKr after Dofetilide block"
            plt.plot(range(len(series)), series, label=self.regimen_label)
        for x in np.arange(0, self.num_hours + 2, tau):
            plt.axvline(x=x, color="r", linestyle="--", linewidth=0.5)
        plt.xlabel("Time (hours)"); plt.ylabel(ylab)
        plt.title(ttl)
        plt.legend(); plt.show()
        self.initial_values()

    # Plot the change of repolzarization duration 
    def plotsimulation_3(self):
        """Simulation 3: t(repolar) (ms)"""
        plt.figure()
        for reg in REGIMENS:
            self.set_regimen(reg["doses"], reg["label"]) ; self.initial_values(); self.run_simulations()
            plt.plot(range(len(self.t_repolar)), self.t_repolar, label=self.regimen_label)
        for x in np.arange(0, self.num_hours + 2, tau):
            plt.axvline(x=x, color="r", linestyle="--", linewidth=0.7)
        plt.xlabel("Time (hours)"); plt.ylabel("t(repolar) (ms)")
        plt.title("Simulation 3: Total AP repolarization time")
        plt.legend(); plt.show()
        self.initial_values()

    # Plot the change of QT interval duration 
    def plotsimulation_4(self):
        """Simulation 4: t(QT) (ms)"""
        plt.figure()
        for reg in REGIMENS:
            self.set_regimen(reg["doses"], reg["label"]) ; self.initial_values(); self.run_simulations()
            plt.plot(range(len(self.t_QT_ms)), self.t_QT_ms, label=self.regimen_label)
        for x in np.arange(0, self.num_hours + 2, tau):
            plt.axvline(x=x, color="r", linestyle="--", linewidth=0.7)
        plt.xlabel("Time (hours)"); plt.ylabel("t(QT) (ms)")
        plt.title("Simulation 4: QT interval duration")
        plt.legend(); plt.show()
        self.initial_values()

    # Helper function to print simulation results in terminal
    def print_parameters(self, summary: bool = True, precision: int = 6):
        """
        Run the simulation and print key parameters for each dosing regimen.

        For each regimen, prints a per-hour table for:
        - C_dof (mg/L), IKr (A), t_repolar (ms), t_QT_ms (ms)
        If summary=True, also prints simple summary stats. State resets per regimen.
        """
        header = (
            f"{'hour':>5}  "
            f"{'C_dof(mg/L)':>14}  "
            f"{'IKr(A)':>14}  "
            f"{'t_repolar(ms)':>14}  "
            f"{'t_QT(ms)':>12}"
        )
        fmt = (
            f"{{hour:5d}}  "
            f"{{C:14.{precision}g}}  "
            f"{{IKr:14.{precision}g}}  "
            f"{{trep:14.{precision}g}}  "
            f"{{qt:12.{precision}g}}"
        )

        for reg in REGIMENS:
            self.set_regimen(reg["doses"], reg["label"]) ; self.initial_values(); self.run_simulations()
            print(f"\n=== {self.regimen_label} ===")
            print(header)
            print('-' * len(header))
            length = min(len(self.C_dof), len(self.IKr), len(self.t_repolar), len(self.t_QT_ms))
            for h in range(length):
                print(fmt.format(hour=h, C=float(self.C_dof[h]), IKr=float(self.IKr[h]), trep=float(self.t_repolar[h]), qt=float(self.t_QT_ms[h])))

            if summary:
                def _stats(arr):
                    arr_np = np.asarray(arr, dtype=float)
                    return float(np.nanmin(arr_np)), float(np.nanmax(arr_np)), float(arr_np[-1])
                C_min, C_max, C_end = _stats(self.C_dof)
                I_min, I_max, I_end = _stats(self.IKr)
                R_min, R_max, R_end = _stats(self.t_repolar)
                Q_min, Q_max, Q_end = _stats(self.t_QT_ms)
                print("Summary (min, max, final):")
                print(f"  C_dof (mg/L):     {C_min:.{precision}g}, {C_max:.{precision}g}, {C_end:.{precision}g}")
                print(f"  IKr (A):          {I_min:.{precision}g}, {I_max:.{precision}g}, {I_end:.{precision}g}")
                print(f"  t_repolar (ms):   {R_min:.{precision}g}, {R_max:.{precision}g}, {R_end:.{precision}g}")
                print(f"  t_QT (ms):        {Q_min:.{precision}g}, {Q_max:.{precision}g}, {Q_end:.{precision}g}")

        self.initial_values()

# Main execution of model
if __name__ == "__main__":
    sim = DosageSimulator(dose_hours=100)
    sim.print_parameters(summary=True, precision=5)
    sim.plotsimulation_1()
    sim.plotsimulation_2()
    sim.plotsimulation_3()
    sim.plotsimulation_4()


dataset_3 = [
             {0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 7, 8, 12, 14, 24},
             {0, 0.842, 1.54, 2.22, 3.19, 3.01, 2.74, 2.21, 2.27, 1.64, 1.75, 1.51, 1.23, 0.694, 0.527, 0.207}
]

