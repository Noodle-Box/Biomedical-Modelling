###################################################################################
# Project File for BIOE3001 Project
# Modelling the effect of Dofetilide on diagnosed Short QT patients
#
# Author: Tevyn Vergara
#
#
# ChatGPT used to restructure in class - based simulation.
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
###################################################################################
### Parameters for Simulation ###

# dosing information
Dosage = [250,125] # dosage of oral Dofetilide (mg) from (Drugs.com, 2025)

# -----------------------------
# Schedule: 250 mg ×2 then 125 mg, q12h
# -----------------------------
tau = 12
t_s = 0
t_e = 2

def dose_at(n: int) -> float:
    return 250 if n < 2 else 125

def delta(t: int, n: int, tau: int) -> int:
    return 1 if (t - n*tau - t_s) == 0 else 0

# -----------------------------
# PK parameters (units in comments)
# -----------------------------
V_SC   = 70.0     # L      (apparent central volume)
KABS   = 0.6     # 1/h    (first-order absorption from gut)
KDISS  = 0.8     # 1      (tablet absorption efficiency)
t_half = 5.0    # h
kDecay  = np.log(2) / t_half  # 1/h

# Extra sink due to binding removing drug from plasma
KBinding = 0.01  # 1/h  (strength of removal; tune to your dataset)

# -----------------------------
# Binding (Eq. 1): Block(C) in [0,1]
# -----------------------------

# Parameters to Tune 
k1 = 2     # Emax (≤1.0), part of Kbinding calculation
k2 = 30     # dimensionless gain on repolarization scaling
k3 = 1.1     # curvature exponent in Eq. 3 (tune)
k4 = 1.1   # scaling from Δt_repolar to QT (ms per ms * factor in Eq. 4)

# Hill Mechanics Parameters
EC50 = 10               # mg/L   (2 µg/L) — tune with your data
hill = 0.015            # Hill n

def hill_mech(v, Emax=k1, EC50=EC50, n=hill):
    v = max(v, 1e-12)
    num = Emax * (v ** n)
    den = (EC50 ** n + v ** n)
    return float(np.clip(num / den, 0.0, 1.0))

# -----------------------------
# Electrophysiology (Eq. 2–4)
# -----------------------------
IKR0_uA       = 20*10**-9  # pA baseline IKr magnitude, from (...)
t0_repolar_ms = 300.0   # ms baseline repolarization time, from (...)
t0_QT_ms      = 350.0   # ms baseline QT, from (...)
EPS = 1e-9

class DosageSimulator:
    """
    Discrete-hour simulator (mirrors your reference structure and the previous code),
    with states and outputs recorded each hour.
    """
    def __init__(self, num_days: int):
        self.num_hours = int(num_days * 24)
        self.initial_values()

    def initial_values(self):

        # Blood drug concentration
        self.C_dof = [0.0]    # mg/L in plasma

        # Potassium current 
        self.IKr = [IKR0_uA]           # µA (post-block current, Eq. 2)

        # Cardiac AP repolarization time
        self.t_repolar = [t0_repolar_ms]

        # QT interval duration
        self.t_QT_ms = [t0_QT_ms]

    def run_simulations(self):
        t = 0
        n = -1  # dose window index

        while t <= self.num_hours:
            if t % tau == 0:
                n += 1 # increment tau to the next dose
            D = dose_at(n)  # mg scheduled for this window (For 250 or 125 mg ingestion)

            """ Simulation 1 - Dofetilide concentration in blood """

            # current plasma level, incoporate heaviside function for dosage timing
            C = float(KDISS*self.C_dof[t])
            C += (D / V_SC) * delta(t, n, tau)

            # over the next hour: elimination + binding sink only (no gut inflow)
            dC = -kDecay * C - KBinding * C
            # advance states (keep A_gut just to preserve structure/length)
            self.C_dof.append(C + dC)


            """ Simulation 2 - IKr after Dofetilide block """

            r_bind = KBinding * self.C_dof[t]
            block = hill_mech(r_bind)          # x in [0,1]
            IKr_current = IKR0_uA * (1.0 - block)
            self.IKr.append(IKr_current)


            """ Simulation 3 - Repolarization time (ms) """

            # --- t_repolar: grows with block; use positive form ---
            #t_repolar = t0_repolar_ms + (k2 * (1 - (IKR0_uA/max(self.IKr[t], EPS))**k3))
            t_repolar = t0_repolar_ms + (k2 * (block))**k3
            self.t_repolar.append(t_repolar)

            
            """ Simulation 4 - QT Interval duration (ms) """
            t_QT = t0_QT_ms + k4 * (self.t_repolar[t] - t0_repolar_ms)
            self.t_QT_ms.append(t_QT)

            t += 1

    # -----------------------------
    # Plotters (your 4 simulations)
    # -----------------------------
    def plot_simulation_1(self):
        """Simulation 1: C(Dof) in blood (mg/L)"""
        self.run_simulations()
        plt.figure()
        plt.plot(range(len(self.C_dof)), self.C_dof, label='250 mg ×2 then 125 mg, q12h')
        for x in np.arange(0, self.num_hours + 2, tau):
            plt.axvline(x=x, color='r', linestyle='--', linewidth=0.7)
        plt.xlabel('Time (hours)'); plt.ylabel('C(Dof) (mg/L)')
        plt.title('Simulation 1 — Dofetilide concentration in blood')
        plt.legend(); plt.show()
        self.initial_values()


    def plot_simulation_2(self, show_reduction=False):
        """Simulation 2: IKr after block (µA). Set show_reduction=True to plot reduction instead."""
        self.run_simulations()
        plt.figure()
        series = self.IKr_reduction_uA if show_reduction else self.IKr
        ylab   = 'IKr reduction (µA)' if show_reduction else 'IKr (µA)'
        plt.plot(range(len(series)), series, label='250 mg ×2 then 125 mg, q12h')
        for x in np.arange(0, self.num_hours + 2, tau):
            plt.axvline(x=x, color='r', linestyle='--', linewidth=0.5)
        plt.xlabel('Time (hours)'); plt.ylabel(ylab)
        ttl = 'Simulation 2 — IKr reduction from Dofetilide' if show_reduction else 'Simulation 2 — IKr after Dofetilide block'
        plt.title(ttl)
        plt.legend(); plt.show()
        self.initial_values()


    def plot_simulation_3(self):
        """Simulation 3: t(repolar) (ms)"""
        self.run_simulations()
        plt.figure()
        plt.plot(range(len(self.t_repolar)), self.t_repolar, label='250 mg ×2 then 125 mg, q12h')
        for x in np.arange(0, self.num_hours + 2, tau):
            plt.axvline(x=x, color='r', linestyle='--', linewidth=0.7)
        plt.xlabel('Time (hours)'); plt.ylabel('t(repolar) (ms)')
        plt.title('Simulation 3 — Total AP repolarization time')
        plt.ylim(bottom=275)
        plt.legend()
        plt.show()
        self.initial_values()


    def plot_simulation_4(self):
        """Simulation 4: t(QT) (ms)"""
        self.run_simulations()
        plt.figure()
        plt.plot(range(len(self.t_QT_ms)), self.t_QT_ms, label='250 mg ×2 then 125 mg, q12h')
        for x in np.arange(0, self.num_hours + 2, tau):
            plt.axvline(x=x, color='r', linestyle='--', linewidth=0.7)
        plt.xlabel('Time (hours)'); plt.ylabel('t(QT) (ms)')
        plt.title('Simulation 4 — QT interval duration')
        plt.ylim(bottom=325)
        plt.legend()
        plt.show()
        self.initial_values()


if __name__ == "__main__":
    sim = DosageSimulator(num_days=4)
    sim.plot_simulation_1()
    sim.plot_simulation_2()
    sim.plot_simulation_3()
    sim.plot_simulation_4()


