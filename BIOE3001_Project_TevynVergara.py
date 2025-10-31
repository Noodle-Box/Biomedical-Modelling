###################################################################################
# Project File for BIOE3001 Project
# Modelling the effect of Dofetilide on diagnosed Short QT patients
#
# Author: Tevyn Vergara
#
# ChatGPT used to restructure in class-based simulation.
###################################################################################
"""
Simulation 1: C(Dof): concentration of Dofetilide in blood (mg/L)
Simulation 2: IKr(Block): Potassium Channel Current (IKr) after block (A)
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
Dosage = [250, 125]  # oral Dofetilide doses (mg)

# Schedule: 250 mg x2 then 125 mg, q12h
tau = 12   # hours between doses
t_s = 0    # dose start offset (h)
t_e = 2    # not used in this simplified model

def dose_at(n: int) -> float:
    return 250.0 if n < 2 else 125.0

def delta(t: int, n: int, tau: int) -> int:
    return 1 if (t - n * tau - t_s) == 0 else 0

# -----------------------------
# PK parameters
# -----------------------------
V_SC   = 70.0     # L, apparent central volume
KABS   = 0.6      # 1/h, first-order absorption from gut (not used directly)
KDISS  = 0.8      # tablet absorption efficiency (dimensionless)
t_half = 5.0      # h
kDecay = np.log(2) / t_half  # 1/h

# Extra sink due to binding removing drug from plasma
KBinding = 0.01  # 1/h

# -----------------------------
# Binding (Eq. 1): Block(C) in [0,1]
# -----------------------------

# Parameters to tune (mechanics from IKr to repolarization/QT)
k1 = 2
k2 = 30      # gain on repolarization scaling
k3 = 1.1     # curvature exponent in Eq. 3
k4 = 1.1     # scaling from t_repolar to QT

# Time‑effective binding parameters (tweakable)\n
# Target rise time to minimum IKr after a dose
# A ramp function scales the effective block to reach its maximum at this time.
T_TO_PEAK_H    = 3.0      # hours to reach IKr minimum (peak block)
RAMP_POWER     = 1.0      # >1 slows early rise; <1 speeds early rise

# Effect-site link
ke0 = 0.15          # 1/h, effect-site rate constant (slower → smaller gradient)
fu_plasma = 0.6     # fraction unbound


""" Hill Function """
# PD Hill on effect-site concentration
Emax_block = 0.9               # max fractional block (leave some residual current)
EC50_C = 0.005                 # mg/L at 50% block
hill_nC = 1.0                  # Hill slope

def hill_on_C(C, Emax=Emax_block, EC50=EC50_C, n=hill_nC):
    C = max(float(C), 1e-12)
    return float(np.clip(Emax * (C**n) / (EC50**n + C**n), 0.0, 1.0))

# -----------------------------
# Electrophysiology (Eq. 2–4)
# -----------------------------
IKR0_uA       = 60e-12     # A baseline IKr magnitude
t0_repolar_ms = 160.0     # ms baseline repolarization time
t0_QT_ms      = 300.0     # ms baseline QT
EPS = 1e-18               # for division safety

class DosageSimulator:
    """
    Discrete-hour simulator with states and outputs recorded each hour.
    """
    def __init__(self, num_days: int):
        self.num_hours = int(num_days * 24)
        self.initial_values()

    def initial_values(self):
        # Blood drug concentration
        self.C_dof = [0.0]    # mg/L in plasma

        # Binding site compartment
        self.Ce = [0.0]       # effect-site conc (mg/L)
        self.block = [0.0]    # fractional block (effective)

        # Potassium current
        self.IKr = [IKR0_uA]  # A (post-block current, Eq. 2)

        # Cardiac AP repolarization time
        self.t_repolar = [t0_repolar_ms]

        # QT interval duration
        self.t_QT_ms = [t0_QT_ms]

    def run_simulations(self):
        t = 0
        n = -1  # dose window index

        while t <= self.num_hours:
            if t % tau == 0:
                n += 1  # new dose window
            D = dose_at(n)  # mg scheduled for this window

            # Simulation 1 — Dofetilide concentration in blood
            C = float(KDISS * self.C_dof[t])
            C += (D / V_SC) * delta(t, n, tau)
            dC = -kDecay * C - KBinding * C
            self.C_dof.append(max(C + dC, 0.0))

            # Simulation 2 — IKr from drug binding (Hill kinetics)
            # Effect-site update
            Ce_prev = self.Ce[t]
            Ce = Ce_prev + ke0 * (C - Ce_prev)  # discrete 1-h Euler
            Ce_u = fu_plasma * Ce
            self.Ce.append(Ce)

            # Fractional block from effect-site concentration
            block = hill_on_C(Ce_u)

            # Start from dose time (no delay/window gating)
            last_dose_h = n * tau + t_s if n >= 0 else -1e9
            start_h = last_dose_h
            block_eff = block

            # Additional ramp to ensure ~T_TO_PEAK_H hours to minimum IKr
            t_center = t + 0.5
            ramp_raw = (t_center - start_h) / max(T_TO_PEAK_H, 1e-9)
            ramp = float(np.clip(ramp_raw, 0.0, 1.0)) ** RAMP_POWER

            block_eff = block_eff * ramp

            IKr_current = IKR0_uA * (1.0 - block_eff)
            self.IKr.append(IKr_current)
            self.block.append(block_eff)

            # Simulation 3 — Repolarization time (ms)
            ratio = IKR0_uA / max(IKr_current, EPS)
            t_repolar = t0_repolar_ms * (1.0 + k2 * (ratio ** k3))
            self.t_repolar.append(t_repolar)

            # Simulation 4 — QT interval duration (ms)
            t_QT = t0_QT_ms + k4 * (t_repolar - t0_repolar_ms)
            self.t_QT_ms.append(t_QT)

            t += 1

    # -----------------------------
    # Plotters
    # -----------------------------
    def plot_simulation_1(self):
        self.run_simulations()
        plt.figure()
        plt.plot(range(len(self.C_dof)), self.C_dof, label='250 mg x2 then 125 mg, q12h')
        for x in np.arange(0, self.num_hours + 2, tau):
            plt.axvline(x=x, color='r', linestyle='--', linewidth=0.7)
        plt.xlabel('Time (hours)'); plt.ylabel('C(Dof) (mg/L)')
        plt.title('Simulation 1 — Dofetilide concentration in blood')
        plt.legend(); plt.show()
        self.initial_values()

    def plot_simulation_2(self, show_reduction=False):
        self.run_simulations()
        plt.figure()
        if show_reduction:
            series = [IKR0_uA - y for y in self.IKr]
            ylab = 'IKr reduction (A)'
            ttl = 'Simulation 2 — IKr reduction from Dofetilide'
        else:
            series = self.IKr
            ylab = 'IKr (A)'
            ttl = 'Simulation 2 — IKr after Dofetilide block'
        plt.plot(range(len(series)), series, label='250 mg x2 then 125 mg, q12h')
        for x in np.arange(0, self.num_hours + 2, tau):
            plt.axvline(x=x, color='r', linestyle='--', linewidth=0.5)
        plt.xlabel('Time (hours)'); plt.ylabel(ylab)
        plt.ylim(bottom = 50e-12 )
        plt.title(ttl)
        plt.legend(); plt.show()
        self.initial_values()

    def plot_simulation_3(self):
        self.run_simulations()
        plt.figure()
        plt.plot(range(len(self.t_repolar)), self.t_repolar, label='250 mg x2 then 125 mg, q12h')
        for x in np.arange(0, self.num_hours + 2, tau):
            plt.axvline(x=x, color='r', linestyle='--', linewidth=0.7)
        plt.xlabel('Time (hours)'); plt.ylabel('t(repolar) (ms)')
        plt.title('Simulation 3 — Total AP repolarization time')
        plt.ylim(bottom=275)
        plt.legend(); plt.show()
        self.initial_values()

    def plot_simulation_4(self):
        self.run_simulations()
        plt.figure()
        plt.plot(range(len(self.t_QT_ms)), self.t_QT_ms, label='250 mg x2 then 125 mg, q12h')
        for x in np.arange(0, self.num_hours + 2, tau):
            plt.axvline(x=x, color='r', linestyle='--', linewidth=0.7)
        plt.xlabel('Time (hours)'); plt.ylabel('t(QT) (ms)')
        plt.title('Simulation 4 — QT interval duration')
        plt.ylim(bottom=325)
        plt.legend(); plt.show()
        self.initial_values()

    def print_parameters(self, summary: bool = True, precision: int = 6):
        """
        Run the simulation and print key parameters to the console.

        Prints a per-hour table for:
          - C_dof (mg/L)
          - IKr (A)
          - t_repolar (ms)
          - t_QT_ms (ms)

        If summary=True, also prints simple summary stats.
        Resets internal state afterward to avoid side-effects on plotting.
        """
        self.run_simulations()
        header = (
            f"{'hour':>5}  "
            f"{'C_dof(mg/L)':>14}  "
            f"{'IKr(A)':>14}  "
            f"{'t_repolar(ms)':>14}  "
            f"{'t_QT(ms)':>12}"
        )
        print(header)
        print('-' * len(header))
        fmt = (
            f"{{hour:5d}}  "
            f"{{C:14.{precision}g}}  "
            f"{{IKr:14.{precision}g}}  "
            f"{{trep:14.{precision}g}}  "
            f"{{qt:12.{precision}g}}"
        )
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
            print("\nSummary (min, max, final):")
            print(f"  C_dof (mg/L):     {C_min:.{precision}g}, {C_max:.{precision}g}, {C_end:.{precision}g}")
            print(f"  IKr (A):          {I_min:.{precision}g}, {I_max:.{precision}g}, {I_end:.{precision}g}")
            print(f"  t_repolar (ms):   {R_min:.{precision}g}, {R_max:.{precision}g}, {R_end:.{precision}g}")
            print(f"  t_QT (ms):        {Q_min:.{precision}g}, {Q_max:.{precision}g}, {Q_end:.{precision}g}")

        self.initial_values()

if __name__ == "__main__":
    sim = DosageSimulator(num_days=4)
    sim.print_parameters()
    sim.plot_simulation_1()
    sim.plot_simulation_2()
    sim.plot_simulation_3()
    sim.plot_simulation_4()

