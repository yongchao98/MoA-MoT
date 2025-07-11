import math

def solve_contaminant_transport():
    """
    Calculates the concentration of three chemicals reaching a spring and assesses the mixture effect.
    """
    # Step 1: Define shared parameters and constants from the problem description
    x_distance = 100  # m
    v_velocity = 0.5  # m/d
    alpha_dispersivity = 0.5  # m
    
    injected_water_volume_L = 1500000  # L
    fracture_length = 80  # m
    fracture_height = 10  # m
    
    coal_density_kg_m3 = 1346  # kg/m³
    foc_fraction_organic_carbon = 0.50  # dimensionless
    n_porosity = 0.30  # dimensionless (from 30% water content)
    
    # Step 2: Calculate derived shared parameters
    # Calculate travel time for water to reach the spring
    t_time = x_distance / v_velocity  # days
    
    # Calculate Dispersion coefficient (D)
    D_dispersion_coeff = alpha_dispersivity * v_velocity  # m²/d
    
    # Calculate bulk density of the coal seam
    rho_b_bulk_density = coal_density_kg_m3 * (1 - n_porosity)  # kg/m³
    
    # Pre-calculate the denominator term from the advection-dispersion equation
    denominator_sqrt = math.sqrt(4 * math.pi * D_dispersion_coeff * t_time)

    # Define chemical properties in a list of dictionaries
    chemicals = [
        {
            'name': 'Atrazine',
            'product_percent': 0.01,
            'conc_in_product_ug_L': 40,
            'log_koc': 2.20,
            'half_life_days': 90,
            'ec50_ug_L': 100
        },
        {
            'name': 'PFOS',
            'product_percent': 0.001,
            'conc_in_product_ug_L': 300,
            'log_koc': 3.65,
            'half_life_days': 14965,
            'ec50_ug_L': 480
        },
        {
            'name': 'Endosulfan',
            'product_percent': 0.005,
            'conc_in_product_ug_L': 20,
            'log_koc': 4.3,
            'half_life_days': 60,
            'ec50_ug_L': 560
        }
    ]

    print("--- Environmental Parameters ---")
    print(f"Distance to spring (x): {x_distance} m")
    print(f"Groundwater velocity (v): {v_velocity} m/d")
    print(f"Travel time (t = x/v): {t_time} days")
    print(f"Dispersion Coefficient (D): {D_dispersion_coeff} m²/d")
    print(f"Coal Bulk Density (ρ_b): {rho_b_bulk_density:.2f} kg/m³")
    print(f"Porosity (n): {n_porosity}")
    print(f"Fraction Organic Carbon (foc): {foc_fraction_organic_carbon}")
    print("-" * 30 + "\n")

    results = []
    total_hazard_index = 0

    # Step 3: Process each chemical
    for chem in chemicals:
        print(f"--- Calculating Concentration for: {chem['name']} ---")
        
        # Calculate total mass injected (M_total)
        m_total_ug = injected_water_volume_L * chem['product_percent'] * chem['conc_in_product_ug_L']
        
        # Calculate mass per unit area (M_area)
        source_area_m2 = fracture_length * fracture_height
        m_area_ug_m2 = m_total_ug / source_area_m2
        
        # Calculate decay rate constant (k)
        k_decay_rate = math.log(2) / chem['half_life_days']
        
        # Calculate Retardation Factor (R)
        koc = 10**chem['log_koc']
        kd_partition_coeff = koc * foc_fraction_organic_carbon # L/kg
        # R = 1 + (ρ_b/n) * Kd. Unit conversion from m³ to L is needed (1m³ = 1000L).
        R_retardation = 1 + (rho_b_bulk_density / n_porosity) * kd_partition_coeff * 0.001

        # The term exp(−(x−vt)^2/4Dt) becomes exp(0) = 1 since we evaluate at t=x/v
        advection_dispersion_term = m_area_ug_m2 / denominator_sqrt
        
        # The term exp(-kt)
        decay_term = math.exp(-k_decay_rate * t_time)
        
        # Retardation term (1/R)
        retardation_term = 1 / R_retardation
        
        # Calculate final concentration in µg/m³
        conc_ug_m3 = advection_dispersion_term * decay_term * retardation_term
        
        # Convert concentration to µg/L (1 m³ = 1000 L)
        conc_ug_L = conc_ug_m3 / 1000.0
        
        # Calculate Hazard Quotient (HQ)
        hq = conc_ug_L / chem['ec50_ug_L']
        total_hazard_index += hq
        
        chem['concentration_ug_L'] = conc_ug_L
        results.append(chem)
        
        print(f"Total Mass (M_total): {m_total_ug:,.2f} µg")
        print(f"Source Area: {source_area_m2} m²")
        print(f"Mass per Area (M_area): {m_area_ug_m2:.2f} µg/m²")
        print(f"Retardation Factor (R): {R_retardation:.2f}")
        print(f"Decay Rate (k): {k_decay_rate:.6f} d⁻¹")
        
        print("\nFinal Concentration Calculation (C = (M_area / sqrt(4πDt)) * (1/R) * exp(-kt)):")
        print(f"C(µg/m³) = ({m_area_ug_m2:.2f} / {denominator_sqrt:.2f}) * ({retardation_term:.6f}) * ({decay_term:.6f})")
        print(f"C = {conc_ug_m3:.8f} µg/m³")
        print(f"C = {conc_ug_L:.8f} µg/L")
        
        print(f"\nHazard Quotient (HQ = C/EC50): {conc_ug_L:.8f} / {chem['ec50_ug_L']} = {hq:.2e}")
        print("-" * 30 + "\n")

    # Step 4: Determine highest concentration and mixture effect
    highest_conc_chem = max(results, key=lambda x: x['concentration_ug_L'])
    
    print("--- FINAL RESULTS ---")
    print(f"The highest concentration reaching the spring is from {highest_conc_chem['name']} at {highest_conc_chem['concentration_ug_L']:.8f} µg/L.")
    
    print(f"\nMixture Effect Analysis:")
    print(f"The total Hazard Index (HI = ΣHQ) is: {total_hazard_index:.2e}")
    if total_hazard_index < 1:
        print("Since the Hazard Index is much less than 1, the risk of an additive effect from the mixture is considered negligible.")
    else:
        print("Since the Hazard Index is greater than or equal to 1, the additive effect of the mixture poses a potential risk to the algal community.")
    print("The model used assumes an additive effect for the mixture components.")


solve_contaminant_transport()
<<<Atrazine>>>