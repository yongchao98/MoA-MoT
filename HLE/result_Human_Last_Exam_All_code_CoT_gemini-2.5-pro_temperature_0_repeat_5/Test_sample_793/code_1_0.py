import math

def solve_contaminant_transport():
    """
    Calculates contaminant concentrations at a spring and assesses mixture effects
    based on the provided hydrogeological problem.
    """
    # --- Step 1: Define initial parameters ---
    # Environmental parameters
    V_water_total_L = 1_500_000.0
    coal_density_kg_m3 = 1346.0  # Bulk density (rho_b)
    foc = 0.50  # Fraction of organic carbon
    porosity = 0.30  # Water content (theta)
    fracture_length_m = 80.0
    fracture_height_m = 10.0
    distance_m = 100.0  # Distance to spring (x)
    advection_velocity_m_d = 0.5  # (v)
    dispersivity_m = 0.5 # (alpha)

    # Chemical properties
    chemicals = [
        {
            "name": "Atrazine",
            "product_percent": 0.01,
            "conc_in_product_ug_L": 40.0,
            "log_koc": 2.20,
            "half_life_d": 90.0,
            "ec50_ug_L": 100.0,
        },
        {
            "name": "PFOS",
            "product_percent": 0.001,
            "conc_in_product_ug_L": 300.0,
            "log_koc": 3.65,
            "half_life_d": 14965.0,
            "ec50_ug_L": 480.0,
        },
        {
            "name": "Endosulfan",
            "product_percent": 0.005,
            "conc_in_product_ug_L": 20.0,
            "log_koc": 4.3,
            "half_life_d": 60.0,
            "ec50_ug_L": 560.0,
        }
    ]

    # --- Step 2: Calculate shared transport parameters ---
    # Dispersion coefficient (D)
    D_m2_d = advection_velocity_m_d * dispersivity_m
    # Time for plume center to reach the spring (t)
    t_d = distance_m / advection_velocity_m_d
    # Cross-sectional area for the initial mass distribution
    A_m2 = fracture_length_m * fracture_height_m

    print("Shared Parameters:")
    print(f"Dispersion Coefficient (D): {D_m2_d:.2f} m^2/d")
    print(f"Advective Travel Time (t): {t_d:.0f} days")
    print(f"Source Area (A): {A_m2:.0f} m^2\n")

    results = []
    hazard_quotients = []

    # --- Step 3: Loop through each chemical and perform calculations ---
    for chem in chemicals:
        print(f"--- Calculating for {chem['name']} ---")

        # a) Calculate initial mass injected (M_total) and mass per area
        V_product_L = V_water_total_L * chem['product_percent']
        M_total_ug = V_product_L * chem['conc_in_product_ug_L']
        M_per_area_ug_m2 = M_total_ug / A_m2
        
        # b) Calculate chemical-specific parameters
        k_per_d = math.log(2) / chem['half_life_d']
        Koc_L_kg = 10**chem['log_koc']
        Kd_L_kg = Koc_L_kg * foc
        # R = 1 + (rho_b / theta) * Kd. Kd needs to be in m^3/kg.
        R = 1 + (coal_density_kg_m3 / porosity) * (Kd_L_kg / 1000.0)

        # c) Apply the Advection-Diffusion Equation
        # C(ug/m3) = (M_per_area / sqrt(4*pi*D*t)) * exp(-kt) * (1/R)
        numerator = M_per_area_ug_m2 * math.exp(-k_per_d * t_d)
        denominator = math.sqrt(4 * math.pi * D_m2_d * t_d) * R
        C_ug_m3 = numerator / denominator
        C_ug_L = C_ug_m3 / 1000.0
        
        print(f"Concentration at Spring: {C_ug_L:.10f} μg/L\n")

        # d) Calculate Hazard Quotient
        hq = C_ug_L / chem['ec50_ug_L']
        hazard_quotients.append(hq)

        results.append({
            "name": chem['name'],
            "concentration_ug_L": C_ug_L,
            "hq": hq
        })

    # --- Step 4: Determine the highest concentration and mixture effect ---
    highest_conc_chem = max(results, key=lambda x: x['concentration_ug_L'])
    hazard_index = sum(hazard_quotients)

    print("--- Final Results ---")
    print(f"The highest concentration reaching the spring is from {highest_conc_chem['name']}.")
    print(f"Highest Concentration: {highest_conc_chem['concentration_ug_L']:.10f} μg/L")

    print("\nMixture Effect Analysis:")
    hi_calc_str = " + ".join([f"{res['hq']:.6e}" for res in results])
    print(f"Hazard Index (HI) = Sum of (C/EC50) = {hi_calc_str} = {hazard_index:.6e}")
    print("Since HI is much less than 1, the risk to the algal community is low.")
    print("Assuming concentration addition (the standard approach), the combined effect of the mixture is additive.\n")

    # --- Step 5: Output the final equation with numbers as requested ---
    # Find the data for the highest concentration chemical
    chem_data = next(c for c in chemicals if c['name'] == highest_conc_chem['name'])
    
    # Recalculate parameters for printing
    M_total_ug = (V_water_total_L * chem_data['product_percent']) * chem_data['conc_in_product_ug_L']
    M_per_area_ug_m2 = M_total_ug / A_m2
    k_per_d = math.log(2) / chem_data['half_life_d']
    Koc_L_kg = 10**chem_data['log_koc']
    Kd_L_kg = Koc_L_kg * foc
    R = 1 + (coal_density_kg_m3 / porosity) * (Kd_L_kg / 1000.0)
    
    print(f"--- Final Equation for {highest_conc_chem['name']} ---")
    print("C(μg/L) = [ (M_total / A) / sqrt(4 * π * D * t) * exp(-k * t) * (1/R) ] / 1000")
    print(f"C = [ ({M_total_ug:,.0f} μg / {A_m2:.0f} m^2) / sqrt(4 * π * {D_m2_d:.2f} m^2/d * {t_d:.0f} d) * exp(-{k_per_d:.6f} 1/d * {t_d:.0f} d) * (1 / {R:.2f}) ] / 1000")
    print(f"C = {highest_conc_chem['concentration_ug_L']:.10f} μg/L")

    # --- Final Answer Block ---
    final_answer_content = f"{highest_conc_chem['concentration_ug_L']:.10f}, additive"
    print(f"\n<<<{final_answer_content}>>>")

solve_contaminant_transport()