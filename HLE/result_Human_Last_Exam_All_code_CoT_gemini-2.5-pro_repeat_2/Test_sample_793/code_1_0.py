import numpy as np
import math

def solve_transport_problem():
    """
    This script calculates contaminant transport and risk based on the provided data.
    """
    # --- Step 1: Define constants and parameters ---
    print("This script calculates the highest concentration of a chemical reaching a spring and assesses the mixture's effect.")
    
    # Physical constants
    x = 100.0  # distance to spring (m)
    v = 0.5    # advection velocity (m/d)
    alpha = 0.5 # dispersivity factor (m)
    D = alpha * v # dispersion coefficient (m^2/d)
    rho_s = 1346.0 # coal seam density (kg/m^3)
    n = 0.30     # porosity
    foc = 0.50     # fraction of organic carbon
    rho_b = (1 - n) * rho_s # bulk density (kg/m^3)

    # Source term calculation
    total_water_vol_L = 1500000.0 # Litres
    # Assume cross-sectional area of aquifer is fracture height (10m) * unit width (1m)
    A = 10.0 # m^2

    print("\n--- Model Parameters ---")
    print(f"Distance to spring (x): {x} m")
    print(f"Groundwater velocity (v): {v} m/day")
    print(f"Dispersion coefficient (D): {D} m^2/day")
    
    # Chemical properties
    chemicals = {
        'Atrazine': {
            'product_percent': 0.01,
            'conc_in_product': 40.0, # ug/L
            'log_koc': 2.20,
            'half_life_d': 90.0,
            'ec50': 100.0 # ug/L
        },
        'PFOS': {
            'product_percent': 0.001,
            'conc_in_product': 300.0, # ug/L
            'log_koc': 3.65,
            'half_life_d': 14965.0,
            'ec50': 480.0 # ug/L
        },
        'Endosulfan': {
            'product_percent': 0.005,
            'conc_in_product': 20.0, # ug/L
            'log_koc': 4.3,
            'half_life_d': 60.0,
            'ec50': 560.0 # ug/L
        }
    }

    # --- Step 2: Calculate chemical-specific parameters ---
    for name, params in chemicals.items():
        # Total Mass M_total (ug)
        product_vol_L = total_water_vol_L * params['product_percent']
        M_total_ug = product_vol_L * params['conc_in_product']
        params['M_total_ug'] = M_total_ug
        
        # Effective mass term for formula M_eff (ug/m^2)
        M_eff = M_total_ug / (n * A)
        params['M_eff'] = M_eff

        # Decay rate k (1/d)
        k = np.log(2) / params['half_life_d']
        params['k'] = k
        
        # Partitioning coefficient Kd (L/kg) and Retardation factor R
        Koc = 10**params['log_koc']
        Kd_L_kg = Koc * foc
        Kd_m3_kg = Kd_L_kg * 0.001 # Convert Kd to m^3/kg
        R = 1 + (rho_b / n) * Kd_m3_kg
        params['R'] = R

    # --- Step 3: Define the concentration function based on the provided formula ---
    def calculate_concentration(t, M_eff, D, v, x, k, R):
        if t <= 0:
            return 0
        
        sqrt_term_val = 4 * math.pi * D * t
        
        # C(x,t) = M_eff / sqrt(4*pi*D*t) * exp(-(x-v*t)^2 / (4*D*t)) * exp(-k*t) * (1/R)
        term1 = M_eff / math.sqrt(sqrt_term_val)
        term2 = math.exp(-(x - v * t)**2 / (4 * D * t))
        term3 = math.exp(-k * t)
        term4 = 1 / R
        
        conc_ug_m3 = term1 * term2 * term3 * term4
        conc_ug_L = conc_ug_m3 / 1000.0 # Convert ug/m^3 to ug/L
        return conc_ug_L

    # --- Step 4: Find the maximum concentration for each chemical by simulating over time ---
    max_concentrations = {}
    time_of_max = {}
    # The peak of the main advection-dispersion term occurs around t = x/v = 200 days.
    # We simulate for a wider range to be sure.
    time_points = np.arange(1, 1001, 1) 

    for name, params in chemicals.items():
        max_c = 0
        max_t = 0
        for t in time_points:
            c = calculate_concentration(t, params['M_eff'], D, v, x, params['k'], params['R'])
            if c > max_c:
                max_c = c
                max_t = t
        max_concentrations[name] = max_c
        time_of_max[name] = max_t

    # --- Step 5: Determine the highest concentration and the corresponding chemical ---
    highest_conc_chemical = max(max_concentrations, key=max_concentrations.get)
    highest_conc_value = max_concentrations[highest_conc_chemical]
    time_of_highest_conc = time_of_max[highest_conc_chemical]

    # --- Step 6: Assess the mixture effect at the time of the peak concentration ---
    concentrations_at_peak_time = {}
    for name, params in chemicals.items():
        concentrations_at_peak_time[name] = calculate_concentration(
            time_of_highest_conc, params['M_eff'], D, v, x, params['k'], params['R']
        )

    # Calculate Hazard Index (HI)
    hazard_index = 0
    hqs = {}
    for name, conc in concentrations_at_peak_time.items():
        hq = conc / chemicals[name]['ec50']
        hqs[name] = hq
        hazard_index += hq

    effect = "additive"
    if hazard_index > 1:
        risk_assessment = f"The Hazard Index (HI) is {hazard_index:.4f}, which is greater than 1, suggesting a potential risk to the algal community assuming an additive effect."
    else:
        risk_assessment = f"The Hazard Index (HI) is {hazard_index:.4f}, which is less than 1, suggesting a low risk to the algal community assuming an additive effect."

    # --- Final Output ---
    print("\n--- Results ---")
    print(f"The highest concentration of an individual chemical reaching the spring is {highest_conc_value:.6f} μg/L.")
    print(f"This chemical is {highest_conc_chemical}.")
    print(f"This peak occurs at approximately {time_of_highest_conc} days.")

    print("\n--- Mixture Effect Analysis ---")
    print(f"The effect of the mixture is assumed to be {effect}.")
    print(f"At the time of the peak ({time_of_highest_conc} days), the concentrations and their Hazard Quotients (HQ = Conc / EC50) are:")
    print(f"  Atrazine: {concentrations_at_peak_time['Atrazine']:.6f} μg/L (HQ = {hqs['Atrazine']:.4f})")
    print(f"  PFOS: {concentrations_at_peak_time['PFOS']:.6f} μg/L (HQ = {hqs['PFOS']:.4f})")
    print(f"  Endosulfan: {concentrations_at_peak_time['Endosulfan']:.6f} μg/L (HQ = {hqs['Endosulfan']:.4f})")
    
    print(f"\nThe Hazard Index (HI) is the sum of the individual HQs:")
    print(f"HI = {hqs['Atrazine']:.4f} + {hqs['PFOS']:.4f} + {hqs['Endosulfan']:.4f} = {hazard_index:.4f}")
    print(risk_assessment)

if __name__ == '__main__':
    solve_transport_problem()