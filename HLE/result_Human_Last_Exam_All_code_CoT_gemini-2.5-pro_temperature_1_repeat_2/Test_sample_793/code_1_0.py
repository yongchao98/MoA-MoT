import math

def solve_environmental_problem():
    """
    Calculates the highest chemical concentration reaching a spring and the mixture effect
    based on the provided hydrogeological and chemical data.
    """

    # --- Environmental Parameters ---
    V_water_injected_L = 1_500_000
    rho_b_kg_m3 = 1346  # Bulk density of coal seam
    foc = 0.50          # Fraction of organic carbon
    n = 0.30            # Porosity / water content
    x_m = 100           # Distance to spring
    v_m_d = 0.5         # Advection velocity
    alpha_m = 0.5       # Dispersivity
    D_m2_d = alpha_m * v_m_d  # Dispersion coefficient
    
    # Cross-sectional area of the fracture's leading edge
    fracture_height_m = 10
    fracture_width_m = 0.005
    A_m2 = fracture_height_m * fracture_width_m

    # --- Chemical Parameters ---
    chemicals = {
        "Atrazine": {
            "product_percent": 0.01,
            "conc_in_product_ug_L": 40,
            "log_koc": 2.20,
            "half_life_d": 90,
            "ec50_ug_L": 100,
        },
        "PFOS": {
            "product_percent": 0.001,
            "conc_in_product_ug_L": 300,
            "log_koc": 3.65,
            "half_life_d": 14965,
            "ec50_ug_L": 480,
        },
        "Endosulfan": {
            "product_percent": 0.005,
            "conc_in_product_ug_L": 20,
            "log_koc": 4.3,
            "half_life_d": 60,
            "ec50_ug_L": 560,
        }
    }

    results = {}
    
    # --- Calculations ---
    # Time for non-retarded plume peak to reach the spring
    t_d = x_m / v_m_d
    
    print("--- Calculating Peak Concentrations at the Spring (x=100m) ---\n")

    for name, params in chemicals.items():
        # 1. Total Mass (MTotal) in micrograms
        m_total_ug = V_water_injected_L * params["product_percent"] * params["conc_in_product_ug_L"]

        # 2. Retardation Factor (R)
        koc_L_kg = 10**params["log_koc"]
        kd_L_kg = koc_L_kg * foc
        # Unit conversion: (kg/m^3) * (L/kg) -> L/m^3. Divide by 1000 to get dimensionless ratio.
        r = 1 + (rho_b_kg_m3 / n) * kd_L_kg / 1000

        # 3. Decay Constant (k) in 1/day
        k_per_d = math.log(2) / params["half_life_d"]

        # 4. Peak Concentration (C_peak)
        # Equation: C = (M_total / (A * R)) * (1 / sqrt(4*pi*D*t)) * exp(-k*t)
        # This occurs at t = x/v because we assume retardation doesn't affect travel time in this model.
        # The exponential term exp(-(x-vt)^2 / 4Dt) becomes exp(0) = 1.
        
        sqrt_term = math.sqrt(4 * math.pi * D_m2_d * t_d)
        exp_term = math.exp(-k_per_d * t_d)
        
        # Concentration in ug/m^3
        c_peak_ug_m3 = (m_total_ug / (A_m2 * r)) * (1 / sqrt_term) * exp_term
        
        # Convert to ug/L (1 m^3 = 1000 L)
        c_peak_ug_L = c_peak_ug_m3 / 1000
        
        results[name] = {
            "C_peak_ug_L": c_peak_ug_L,
            "EC50": params["ec50_ug_L"]
        }
        
        # Print the calculation steps
        print(f"Calculating for {name}:")
        print(f"C = (M_total / (A * R)) * (1 / sqrt(4*pi*D*t)) * exp(-k*t)")
        print(f"C = ({m_total_ug:.0f} ug / ({A_m2} m^2 * {r:.2f})) * (1 / sqrt(4*pi*{D_m2_d}*{t_d})) * exp(-{k_per_d:.6f}*{t_d})")
        print(f"C = {c_peak_ug_L:.6f} ug/L\n")

    # --- Find Highest Concentration ---
    highest_conc = 0
    highest_conc_chem = None
    for name, res in results.items():
        if res["C_peak_ug_L"] > highest_conc:
            highest_conc = res["C_peak_ug_L"]
            highest_conc_chem = name
            
    print(f"--- Highest Individual Chemical Concentration ---")
    print(f"The highest concentration is for {highest_conc_chem} with a value of {highest_conc:.6f} ug/L.\n")
    
    # --- Assess Mixture Effect ---
    hqs = {}
    hi = 0
    print("--- Mixture Effect Assessment (Hazard Index) ---")
    print("HI = sum(Concentration / EC50)")
    hi_calc_str = "HI = "
    for name, res in results.items():
        hq = res["C_peak_ug_L"] / res["EC50"]
        hqs[name] = hq
        hi += hq
        hi_calc_str += f"({res['C_peak_ug_L']:.6f}/{res['EC50']}) + "
    
    hi_calc_str = hi_calc_str[:-3] # Remove trailing " + "
    print(hi_calc_str)
    print(f"HI = {hqs['Atrazine']:.6f} (Atrazine) + {hqs['PFOS']:.6f} (PFOS) + {hqs['Endosulfan']:.6f} (Endosulfan)")
    print(f"Total Hazard Index (HI) = {hi:.6f}\n")
    
    print("--- Final Conclusion ---")
    print(f"The highest concentration of an individual chemical reaching the spring is {highest_conc:.6f} μg/L ({highest_conc_chem}).")
    print("The method of summing Hazard Quotients (HQs) to get a Hazard Index (HI) assumes an ADDITIVE effect for the chemical mixture.")


solve_environmental_problem()
<<<The highest concentration is 0.287635 μg/L for Atrazine, and the mixture effect is assumed to be additive.>>>