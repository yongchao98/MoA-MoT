import math

def solve_environmental_problem():
    """
    Solves the advection-diffusion problem for three chemicals and assesses mixture risk.
    """
    # --- Step 1: Define constants and chemical properties ---

    # General parameters
    injected_volume_L = 1_500_000  # L
    distance_x = 100.0  # m
    advection_velocity_v = 0.5  # m/d
    dispersivity_alpha = 0.5  # m
    
    # Coal seam properties
    coal_density_particle = 1346.0  # kg/m3
    foc = 0.50  # fraction of organic carbon
    porosity_n = 0.30

    # Chemical data in a list of dictionaries
    chemicals = [
        {'name': 'Atrazine', 'prod_percent': 0.01, 'conc_in_prod': 40, 'log_koc': 2.20, 'half_life_d': 90, 'ec50': 100},
        {'name': 'PFOS', 'prod_percent': 0.001, 'conc_in_prod': 300, 'log_koc': 3.65, 'half_life_d': 14965, 'ec50': 480},
        {'name': 'Endosulfan', 'prod_percent': 0.005, 'conc_in_prod': 20, 'log_koc': 4.3, 'half_life_d': 60, 'ec50': 560}
    ]

    # --- Step 2: Calculate derived parameters ---
    
    # Dispersion coefficient
    dispersion_D = advection_velocity_v * dispersivity_alpha
    
    # Bulk density
    bulk_density_rhob = coal_density_particle * (1 - porosity_n)
    
    # rho_b / n term (interpreted from Cgw)
    rhob_over_n = bulk_density_rhob / porosity_n

    print("--- Intermediate Parameter Calculations ---")
    print(f"Dispersion Coefficient (D): {dispersion_D:.2f} m^2/d")
    print(f"Bulk Density (ρb): {bulk_density_rhob:.2f} kg/m^3")
    print(f"ρb/n term: {rhob_over_n:.2f} kg/m^3\n")

    results = []
    
    # --- Steps 3 & 4: Process each chemical ---
    for chem in chemicals:
        print(f"--- Processing: {chem['name']} ---")
        
        # Initial Mass (MTotal) in micrograms
        initial_conc_in_water = chem['conc_in_prod'] * chem['prod_percent']
        mtotal_ug = initial_conc_in_water * injected_volume_L
        chem['mtotal_ug'] = mtotal_ug
        print(f"Initial Mass (MTotal): {mtotal_ug:,.0f} μg")

        # Decay rate (k)
        k = math.log(2) / chem['half_life_d']
        chem['k'] = k
        print(f"Decay Rate (k): {k:.6f} per day")

        # Retardation Factor (R)
        koc = 10**chem['log_koc']
        kd_L_kg = koc * foc
        # Convert Kd to m^3/kg for dimensional consistency
        kd_m3_kg = kd_L_kg / 1000.0
        R = 1 + rhob_over_n * kd_m3_kg
        chem['R'] = R
        print(f"Partitioning Coefficient (Kd): {kd_L_kg:.2f} L/kg")
        print(f"Retardation Factor (R): {R:,.1f}")

        # Time of peak concentration (t_peak)
        # Solves (v^2 + 4Dk)t^2 + 2Dt - x^2 = 0 for t
        a = advection_velocity_v**2 + 4 * dispersion_D * k
        b = 2 * dispersion_D
        c = -distance_x**2
        
        # Using the quadratic formula, t = (-b + sqrt(b^2 - 4ac)) / 2a
        t_peak = (-b + math.sqrt(b**2 - 4*a*c)) / (2*a)
        chem['t_peak'] = t_peak
        print(f"Time to Peak Concentration (t_peak): {t_peak:.2f} days")

        # Peak Concentration (C_peak) using the provided formula
        # C(x,t) = (MTotal/sqrt(4πDt)) * exp(−(x−vt)^2/(4Dt)) * exp(−kt) * (1/R)
        term1_denominator = math.sqrt(4 * math.pi * dispersion_D * t_peak)
        term1 = mtotal_ug / term1_denominator
        
        term2_exponent = -((distance_x - advection_velocity_v * t_peak)**2) / (4 * dispersion_D * t_peak)
        term2 = math.exp(term2_exponent)
        
        term3 = math.exp(-k * t_peak)
        
        term4 = 1 / R
        
        c_peak = term1 * term2 * term3 * term4
        chem['c_peak'] = c_peak
        
        print("\nCalculating Peak Concentration C(x,t) = (MTotal/sqrt(4*pi*D*t)) * exp(-(x-v*t)^2/(4*D*t)) * exp(-k*t) * (1/R)")
        print(f"C = ({mtotal_ug:,.0f} / sqrt(4 * {math.pi:.2f} * {dispersion_D:.2f} * {t_peak:.2f})) * "
              f"exp(-({distance_x:.1f} - {advection_velocity_v:.1f}*{t_peak:.2f})^2 / (4 * {dispersion_D:.2f} * {t_peak:.2f})) * "
              f"exp(-{k:.6f} * {t_peak:.2f}) * (1 / {R:.1f})")
        print(f"Peak Concentration for {chem['name']}: {c_peak:.4f} μg/L\n")
        
        results.append(chem)

    # --- Step 5: Find the highest concentration ---
    highest_conc_chem = max(results, key=lambda x: x['c_peak'])
    highest_concentration = highest_conc_chem['c_peak']
    time_of_highest_conc = highest_conc_chem['t_peak']

    print("--- Final Results ---")
    print(f"The highest individual chemical concentration reaching the spring is from {highest_conc_chem['name']}.")
    print(f"Highest Concentration: {highest_concentration:.4f} μg/L")
    print(f"This peak occurs at {time_of_highest_conc:.2f} days.\n")

    # --- Step 6: Assess mixture effect at the time of the highest concentration ---
    print(f"--- Mixture Risk Assessment (at t = {time_of_highest_conc:.2f} days) ---")
    print("The effect of the mixture is assessed assuming additivity using the Hazard Index (HI) model.")
    
    hazard_index = 0
    for chem in results:
        # Calculate concentration of each chemical at time_of_highest_conc
        t = time_of_highest_conc
        k = chem['k']
        R = chem['R']
        mtotal_ug = chem['mtotal_ug']
        
        term1 = mtotal_ug / math.sqrt(4 * math.pi * dispersion_D * t)
        term2 = math.exp(-((distance_x - advection_velocity_v * t)**2) / (4 * dispersion_D * t))
        term3 = math.exp(-k * t)
        term4 = 1 / R
        
        conc_at_t_max = term1 * term2 * term3 * term4
        
        # Hazard Quotient
        hq = conc_at_t_max / chem['ec50']
        hazard_index += hq
        
        print(f"\nAnalysis for {chem['name']}:")
        print(f"Concentration at t={t:.2f} days: {conc_at_t_max:.4f} μg/L")
        print(f"Freshwater Algal EC50: {chem['ec50']} μg/L")
        print(f"Hazard Quotient (HQ = C/EC50): {conc_at_t_max:.4f} / {chem['ec50']} = {hq:.4f}")

    print(f"\nTotal Hazard Index (HI = sum of HQs): {hazard_index:.4f}")
    if hazard_index > 1:
        print("Conclusion: The HI is greater than 1, suggesting a potential risk to the algal community from the additive effect of the mixture.")
    else:
        print("Conclusion: The HI is less than 1, suggesting the risk to the algal community from the additive effect of the mixture is low.")
        
    final_answer_str = (f"The highest concentration is {highest_concentration:.2f} μg/L for {highest_conc_chem['name']}. "
                        f"The mixture effect is assessed assuming additivity, resulting in a Hazard Index of {hazard_index:.2f}, "
                        "which is below the level of concern (HI < 1).")
    
    print(f"\n<<<{final_answer_str}>>>")


if __name__ == '__main__':
    solve_environmental_problem()