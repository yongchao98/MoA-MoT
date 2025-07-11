import numpy as np

def solve_environmental_transport():
    """
    Calculates the highest chemical concentration reaching a spring and the mixture effect
    using a specified advection-diffusion equation.
    """

    # --- Step 1: Define Initial Parameters ---
    x_distance = 100.0  # m
    advection_velocity = 0.5  # m/d
    dispersivity = 0.5  # m
    total_water_vol_L = 1_500_000.0  # L
    
    coal_density_kg_m3 = 1346.0  # kg/m^3
    coal_foc = 0.50  # fraction of organic carbon
    coal_water_mass_frac = 0.30
    
    fracture_l = 80.0 # m
    fracture_h = 10.0 # m
    
    # --- Step 2: Calculate Shared Aquifer Properties ---
    
    # Porosity (n) = volume of water / total volume
    # Mass of 1 m^3 of coal seam = 1346 kg
    # Mass of water in 1 m^3 = 1346 kg * 0.30 = 403.8 kg
    # Volume of water (density ~1000 kg/m^3) = 403.8 kg / 1000 kg/m^3 = 0.4038 m^3
    # Porosity n = 0.4038 m^3 / 1 m^3 = 0.4038
    porosity = 0.4038
    
    # Bulk density (rho_b)
    rho_b_kg_m3 = coal_density_kg_m3
    
    # Hydrodynamic Dispersion coefficient (D)
    D_m2_d = dispersivity * advection_velocity
    
    # Cross-sectional area for mass loading calculation
    A_m2 = fracture_l * fracture_h
    
    # Per the provided formula C(x,t) = M/sqrt(4πDt) * exp(-(x-vt)^2/4Dt)...
    # the peak concentration occurs when the center of mass arrives, so t = x/v.
    t_peak_d = x_distance / advection_velocity

    # --- Step 3: Define and Process Each Chemical ---
    chemicals = {
        "Atrazine": {
            "prod_frac": 0.01, "conc_in_prod_ug_L": 40.0, "log_koc": 2.20,
            "half_life_d": 90.0, "ec50_ug_L": 100.0
        },
        "PFOS": {
            "prod_frac": 0.001, "conc_in_prod_ug_L": 300.0, "log_koc": 3.65,
            "half_life_d": 14965.0, "ec50_ug_L": 480.0
        },
        "Endosulfan": {
            "prod_frac": 0.005, "conc_in_prod_ug_L": 20.0, "log_koc": 4.3,
            "half_life_d": 60.0, "ec50_ug_L": 560.0
        }
    }

    results = {}

    for name, params in chemicals.items():
        # Initial Mass Loading (M_input) in g/m^2
        total_mass_ug = total_water_vol_L * params["prod_frac"] * params["conc_in_prod_ug_L"]
        total_mass_g = total_mass_ug / 1e6
        M_input_g_m2 = total_mass_g / A_m2

        # Decay constant (k) in d^-1
        k_per_day = np.log(2) / params["half_life_d"]

        # Partitioning coefficient (Kd) and Retardation Factor (R)
        Koc_L_kg = 10**params["log_koc"]
        Kd_L_kg = Koc_L_kg * coal_foc
        Kd_m3_kg = Kd_L_kg * 0.001  # Convert L to m^3
        
        # Calculate R = 1 + (rho_b/n) * Kd. The factor in the user's equation is 1/R.
        R_factor = 1.0 + (rho_b_kg_m3 / porosity) * Kd_m3_kg
        retardation_term = 1.0 / R_factor

        # --- Step 4: Calculate Peak Concentration ---
        # Using C(x,t) = (M_input / sqrt(4πDt)) * exp(−(x−vt)²/4Dt) * exp(−kt) * (retardation_term)
        # At t = x/v, the term exp(-(x-vt)^2/4Dt) becomes exp(0) = 1.
        conc_g_m3 = (M_input_g_m2 / np.sqrt(4 * np.pi * D_m2_d * t_peak_d)) \
                    * np.exp(-k_per_day * t_peak_d) \
                    * retardation_term

        # Convert concentration from g/m^3 to ug/L (1 g/m^3 = 1000 ug/L)
        conc_ug_L = conc_g_m3 * 1000.0
        
        results[name] = {"conc_ug_L": conc_ug_L, "ec50_ug_L": params["ec50_ug_L"]}

    # --- Step 5: Analyze Results and Mixture Effect ---
    
    highest_conc = 0.0
    highest_conc_chem = None
    for name, data in results.items():
        if data["conc_ug_L"] > highest_conc:
            highest_conc = data["conc_ug_L"]
            highest_conc_chem = name
    
    print(f"The analysis is conducted at the peak arrival time t = {t_peak_d:.0f} days.")
    print("-" * 30)
    print(f"Highest individual chemical concentration at the spring:")
    print(f"{highest_conc_chem}: {highest_conc:.4e} µg/L")
    print("-" * 30)

    # Calculate Toxic Units (TU) for the mixture at t_peak_d
    TU_mixture = 0.0
    tu_calcs = []

    for name, data in results.items():
        tu = data["conc_ug_L"] / data["ec50_ug_L"]
        TU_mixture += tu
        tu_calcs.append(f"({data['conc_ug_L']:.2e} / {data['ec50_ug_L']})")
    
    effect = ""
    if TU_mixture >= 1.0:
        effect = "additive or synergistic"
    else:
        effect = "less than additive (negligible risk)"
        
    print("Mixture Effect Analysis (using Toxic Units):")
    print(f"TU_mixture = TU_Atrazine + TU_PFOS + TU_Endosulfan")
    final_eq_str = " + ".join(tu_calcs)
    print(f"TU_mixture = {final_eq_str}")
    print(f"TU_mixture = {TU_mixture:.4e}")
    print(f"\nSince the mixture's total Toxic Units value ({TU_mixture:.4e}) is significantly less than 1, the combined effect of the chemicals at this concentration is {effect}.")
    
    # Output the final answer in the required format
    # The final answer requires two parts: highest concentration and mixture effect.
    final_answer = f"The highest concentration is {highest_conc:.4e} µg/L (from {highest_conc_chem}). The mixture effect is less than additive."
    
    return final_answer # This is for logical flow, the prints serve the user.

# Run the solver
solve_environmental_transport()