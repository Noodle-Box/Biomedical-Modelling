## Sample Calculator ##

MM_Dofetillide = 44.157  # g/mol

Dosing_Regimens = {
    0: {"desc": "First Dose: 500mcg then 250mcg", "doses_mcg": [500, 250]},
    1: {"desc": "First Dose: 250mcg then 125mcg", "doses_mcg": [250, 125]},
    2: {"desc": "First Dose: 125mcg then 125mcg", "doses_mcg": [125, 125]},
}

def concentrations_nM(molar_mass_g_per_mol, regimens, volume_L=70.0):
    """
    Compute concentrations (nM) for each dose in regimens.
    - molar_mass_g_per_mol: g/mol
    - regimens: dict with 'doses_mcg' lists (micrograms)
    - volume_L: volume in litres to assume for dilution (default 1 L)
    """
    results = {}
    for key, info in regimens.items():
        concs = []
        for dose_mcg in info.get("doses_mcg", []):
            dose_g = dose_mcg * 1e-9                  # mcg -> g
            molarity_M = dose_g / molar_mass_g_per_mol / volume_L
            conc_nM = molarity_M * 1e9                # M -> nM
            concs.append(conc_nM)
        results[key] = {"desc": info.get("desc",""), "concentrations_nM": concs}
    return results

# Example usage (specify a realistic volume, e.g. 3 L for central compartment):
if __name__ == "__main__":
    res = concentrations_nM(MM_Dofetillide, Dosing_Regimens, volume_L=3.0)
    for k, v in res.items():
        print(k, v["desc"], "->", ["{:.2f} uM".format(x) for x in v["concentrations_nM"]])