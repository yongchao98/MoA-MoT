import math
import numpy as np

def solve_contaminant_transport():
    """
    This script calculates contaminant transport from a hydraulic fracturing site
    to a nearby spring using an advection-diffusion model. It finds the peak
    concentration for three chemicals (Atrazine, PFOS, Endosulfan) and assesses
    the combined mixture effect at that time.
    """
    
    # --- 1. Define constants and initial values from the problem statement ---
    V_water_L = 1500000.0  # Total injected water in litres
    rho_b_kg_m3 = 1346.0   # Coal seam density in kg/m^3 (assumed as bulk density)
    foc = 0.50             # Fraction of organic carbon in coal
    n = 0.30               # Porosity of the coal seam (from water content)
    x_m = 100.0            # Distance to the spring in meters
    v_m_d = 0.5            # Advection velocity in meters/day
    alpha_L_m = 0.5        # Dispersivity factor in meters
    
    # Fracture dimensions used to estimate flow area
    frac_h_m = 10.0
    frac_w_m = 0.005

    # --- 2. Calculate non-chemical-specific derived parameters ---
    # Dispersion coefficient (D) in m^2/day
    D_m2_d = alpha_L_m * v_m_d
    # Cross-sectional area (A) of flow, assuming it's the fracture face
    A_m2 = frac_h_m * frac_w_m
    # Density/porosity term for retardation, converting m^3 to L for unit consistency
    rho_b_n_converted_kg_L = (rho_b_kg_m3 / n) * (1.0 / 1000.0)

    # --- 3. Define chemical-specific data ---
    chemicals = {
        'Atrazine': {
            'product_percent': 0.01, 'C_init_ug_L': 40.0, 'logKoc': 2.20,
            'half_life_d': 90.0, 'EC50_ug_L': 100.0,
        },
        'PFOS': {
            'product_percent': 0.001, 'C_init_ug_L': 300.0, 'logKoc': 3.65,
            'half_life_d': 14965.0, 'EC50_ug_L': 480.0,
        },
        'Endosulfan': {
            'product_percent': 0.005, 'C_init_ug_L': 20.0, 'logKoc': 4.30,
            'half_life_d': 60.0, 'EC50_ug_L': 560.0,
        }
    }

    # --- 4. Calculate derived parameters for each chemical ---
    for name, params in chemicals.items():
        params['MTotal_ug'] = V_water_L * params['product_percent'] * params['C_init_ug_L']
        params['Koc_L_kg'] = 10**params['logKoc']
        params['Kd_L_kg'] = params['Koc_L_kg'] * foc
        params['R'] = 1.0 + rho_b_n_converted_kg_L * params['Kd_L_kg']
        # Add a check for zero half-life to avoid division by zero
        if params['half_life_d'] > 0:
            params['k_d'] = math.log(2) / params['half_life_d']
        else:
            params['k_d'] = 0

    # --- 5. Define the concentration function based on a standard 1D equation ---
    def calculate_concentration_ug_L(t, params):
        if t <= 0: return 0.0
        
        MTotal = params['MTotal_ug']
        R = params['R']
        k = params['k_d']
        
        sqrt_term = math.sqrt(4 * math.pi * D_m2_d * t)
        if sqrt_term == 0: return 0.0
        pre_exp_factor = MTotal / (A_m2 * n * sqrt_term)

        exp_term1_numerator = -((x_m - v_m_d * t)**2)
        exp_term1_denominator = 4 * D_m2_d * t
        exp_term1 = math.exp(exp_term1_numerator / exp_term1_denominator)

        exp_term2 = math.exp(-k * t)
        retardation_factor_term = 1.0 / R

        # Concentration in ug/m^3, then converted to ug/L
        concentration_ug_m3 = pre_exp_factor * exp_term1 * exp_term2 * retardation_factor_term
        return concentration_ug_m3 / 1000.0

    # --- 6. Simulate over time to find peak concentrations ---
    time_days = np.arange(1, 1001, 1)  # Simulate for 1000 days
    peak_concentrations = {}
    peak_times = {}

    for name, params in chemicals.items():
        concentrations = [calculate_concentration_ug_L(t, params) for t in time_days]
        max_c = max(concentrations)
        max_t = time_days[np.argmax(concentrations)]
        peak_concentrations[name] = max_c
        peak_times[name] = max_t

    # --- 7. Identify the highest overall concentration ---
    highest_conc = 0
    highest_conc_chemical = None
    for name, conc in peak_concentrations.items():
        if conc > highest_conc:
            highest_conc = conc
            highest_conc_chemical = name

    # --- 8. Calculate mixture effect at the time of the highest peak ---
    t_of_highest_peak = peak_times[highest_conc_chemical]
    concentrations_at_peak_time = {name: calculate_concentration_ug_L(t_of_highest_peak, params)
                                   for name, params in chemicals.items()}
    hazard_index = sum(c / chemicals[name]['EC50_ug_L'] for name, c in concentrations_at_peak_time.items())
    
    if hazard_index > 1:
        mixture_effect_desc = f"additive, and poses a potential risk (Hazard Index = {hazard_index:.4f})"
    else:
        mixture_effect_desc = f"additive, but the combined risk is likely low (Hazard Index = {hazard_index:.4f})"
    
    final_effect_type = "additive"

    # --- Final Output ---
    print("--- Contaminant Transport Analysis ---")

    # Display the numbers for the equation of the highest concentration chemical
    print(f"\nCalculation breakdown for the highest peak ({highest_conc_chemical}):")
    p = chemicals[highest_conc_chemical]
    t = peak_times[highest_conc_chemical]
    print(f"The peak concentration is calculated at time t = {t} days.")
    print(f"C(t={t}) = (MTotal / (A * n * sqrt(4*pi*D*t))) * exp(-(x-v*t)^2/(4*D*t)) * exp(-k*t) * (1/R)")
    print(f"C({t}) = ({p['MTotal_ug']:.2e} / ({A_m2} * {n} * sqrt(4*pi*{D_m2_d}*{t}))) * exp(-({x_m}-{v_m_d}*{t})^2/(4*{D_m2_d}*{t})) * exp(-{p['k_d']:.4f}*{t}) * (1/{p['R']:.2f})")
    print(f"Resulting Peak Concentration: {highest_conc:.4f} μg/L")
    
    print("\n--- Results Summary ---")
    print(f"The highest concentration reached by an individual chemical is {highest_conc:.4f} μg/L, from {highest_conc_chemical}.")
    print(f"Based on the Concentration Addition model, the predicted mixture effect is {mixture_effect_desc}.")

    final_answer = f"The highest concentration of an individual chemical reaching the spring is {highest_conc:.4f} μg/L. The effect of the mixture at this concentration is {final_effect_type}."
    print(f"\n{final_answer}")
    print(f"\n<<<{highest_conc:.4f}, {final_effect_type}>>>")

if __name__ == '__main__':
    solve_contaminant_transport()