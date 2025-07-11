import numpy as np

def solve_contamination_scenario():
    """
    Calculates the highest chemical concentration reaching a spring and the mixture effect
    based on the Advection-Diffusion Equation and provided site parameters.
    """

    # --- 1. Define Initial Parameters and Constants ---
    # Environmental parameters
    params = {
        'total_water_injected_L': 1500000,
        'coal_density_kg_m3': 1346,
        'organic_carbon_frac': 0.50,
        'water_content_frac': 0.30,  # Porosity (n)
        'fracture_h_m': 10,
        'fracture_w_m': 0.005,
        'distance_to_spring_m': 100,
        'advection_velocity_m_d': 0.5,
        'dispersivity_m': 0.5,
    }

    # Chemical-specific parameters
    chemicals = {
        'Atrazine': {
            'product_percent': 0.01,
            'conc_in_product_ug_L': 40,
            'log_koc': 2.20,
            'half_life_d': 90,
            'ec50_ug_L': 100,
        },
        'PFOS': {
            'product_percent': 0.001,
            'conc_in_product_ug_L': 300,
            'log_koc': 3.65,
            'half_life_d': 14965,
            'ec50_ug_L': 480,
        },
        'Endosulfan': {
            'product_percent': 0.005,
            'conc_in_product_ug_L': 20,
            'log_koc': 4.30,
            'half_life_d': 60,
            'ec50_ug_L': 560,
        }
    }

    # --- 2. Calculate Derived Parameters ---
    
    # Convert total injected water from Litres to m^3
    total_water_m3 = params['total_water_injected_L'] / 1000

    # Calculate dispersion coefficient (D)
    D = params['dispersivity_m'] * params['advection_velocity_m_d']
    
    # Calculate coal seam bulk density (rhob), used for Cgw in the formula
    # rhob = solid_density * (1 - porosity)
    rhob_kg_m3 = params['coal_density_kg_m3'] * (1 - params['water_content_frac'])
    
    # Calculate cross-sectional area (A) of the fracture perpendicular to flow
    A_m2 = params['fracture_h_m'] * params['fracture_w_m']
    
    # Porosity (n)
    n = params['water_content_frac']

    # Process each chemical to calculate its specific equation parameters
    for name, prop in chemicals.items():
        # Total mass injected in kg
        # M_actual = (Total_Water * Product %) * Concentration
        m_actual_kg = (total_water_m3 * prop['product_percent']) * (prop['conc_in_product_ug_L'] * 1e-9)
        prop['m_actual_kg'] = m_actual_kg

        # Normalized mass loading (MTotal) for the equation to make it dimensionally correct
        # MTotal = M_actual / (A * n)
        prop['MTotal_kg_m2'] = m_actual_kg / (A_m2 * n)
        
        # First-order decay rate (k)
        prop['k_per_day'] = np.log(2) / prop['half_life_d']
        
        # Partitioning coefficient (Kd)
        # Koc is in L/kg. Convert to m^3/kg (1 L = 0.001 m^3)
        koc_m3_kg = (10**prop['log_koc']) / 1000
        prop['kd_m3_kg'] = koc_m3_kg * params['organic_carbon_frac']

    # --- 3. Find Peak Concentration for each chemical ---

    # Define the Advection-Diffusion Equation from the prompt
    # Note: The formula's retardation term is a simple multiplier, which is physically non-standard
    # but implemented as per the prompt.
    def get_concentration_kg_m3(t, chem_props):
        MTotal = chem_props['MTotal_kg_m2']
        k = chem_props['k_per_day']
        Kd = chem_props['kd_m3_kg']
        x = params['distance_to_spring_m']
        v = params['advection_velocity_m_d']
        Cgw = rhob_kg_m3 # Using bulk density for Cgw

        # Avoid division by zero at t=0
        if t == 0:
            return 0
        
        sqrt_term = np.sqrt(4 * np.pi * D * t)
        
        # Handle potential overflow in exp by checking the argument
        exp_arg = -((x - v * t)**2) / (4 * D * t)
        if exp_arg < -700: # Below approx exp(-700), the result is effectively zero
            return 0

        transport_term = (MTotal / sqrt_term) * np.exp(exp_arg)
        decay_term = np.exp(-k * t)
        retardation_term = 1 / (1 + Kd * Cgw)
        
        return transport_term * decay_term * retardation_term

    # Simulate over time to find the peak
    # Peak should be around t = x/v = 100/0.5 = 200 days. We'll check up to 1000 days.
    time_days = np.linspace(0.1, 1000, 10000)
    
    peak_info = {}
    for name, prop in chemicals.items():
        concentrations_kg_m3 = np.array([get_concentration_kg_m3(t, prop) for t in time_days])
        
        # Convert kg/m^3 to ug/L for final output (1 kg/m^3 = 1,000,000 ug/L)
        concentrations_ug_L = concentrations_kg_m3 * 1e6

        peak_conc_ug_L = np.max(concentrations_ug_L)
        peak_time_d = time_days[np.argmax(concentrations_ug_L)]
        
        peak_info[name] = {
            'peak_conc_ug_L': peak_conc_ug_L,
            'peak_time_d': peak_time_d,
        }

    # --- 4. Determine Highest Concentration and Mixture Effect ---

    # Find the chemical with the highest peak concentration
    highest_conc_chem = max(peak_info, key=lambda k: peak_info[k]['peak_conc_ug_L'])
    highest_conc_val = peak_info[highest_conc_chem]['peak_conc_ug_L']
    time_of_highest_conc = peak_info[highest_conc_chem]['peak_time_d']

    print(f"Analysis Results:\n")
    print(f"The highest individual concentration is from '{highest_conc_chem}'.")
    print(f"Peak Concentration: {highest_conc_val:.4f} μg/L")
    print(f"This peak occurs at approximately {time_of_highest_conc:.1f} days.\n")
    
    # Calculate the mixture effect at the time the highest concentration occurs
    print(f"Calculating mixture effect at t = {time_of_highest_conc:.1f} days...")
    
    hazard_index = 0
    rq_calcs = []
    for name, prop in chemicals.items():
        # Calculate concentration of each chemical at the specific peak time
        conc_at_peak_time_kg_m3 = get_concentration_kg_m3(time_of_highest_conc, prop)
        conc_at_peak_time_ug_L = conc_at_peak_time_kg_m3 * 1e6
        
        # Calculate Risk Quotient (RQ = C / EC50)
        rq = conc_at_peak_time_ug_L / prop['ec50_ug_L']
        hazard_index += rq
        
        rq_calcs.append(f"  - RQ for {name}: {conc_at_peak_time_ug_L:.4f} μg/L / {prop['ec50_ug_L']} μg/L = {rq:.4f}")

    print("Risk Quotients (RQ) for each chemical:")
    for calc_str in rq_calcs:
        print(calc_str)
        
    print(f"\nHazard Index (HI) = Sum of RQs = {hazard_index:.4f}")

    if hazard_index > 1:
        effect_conclusion = "The mixture effect is additive, and the Hazard Index > 1 indicates a potential risk to the algal community."
    else:
        effect_conclusion = "The mixture effect is additive, and the Hazard Index < 1 suggests the risk to the algal community is likely low."

    print(f"\nConclusion: {effect_conclusion}")

    final_answer = f"The highest concentration of an individual chemical reaching the spring is {highest_conc_val:.4f} μg/L (from {highest_conc_chem}). The effect of the mixture at this concentration is additive, with a calculated Hazard Index of {hazard_index:.4f}, which indicates a potential risk to the algal community."
    print(f"\n<<<{final_answer}>>>")

# Execute the main function
solve_contamination_scenario()