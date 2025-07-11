import math

def solve_environmental_problem():
    """
    Solves the problem by calculating peak chemical concentrations and assessing mixture toxicity.
    """
    
    # --- Step 1: Define constants and parameters ---
    # Environmental parameters
    total_water_volume_L = 1500000.0
    coal_seam_density_kg_m3 = 1346.0
    coal_foc = 0.50  # Fraction of organic carbon
    x_distance_m = 100.0
    v_velocity_m_d = 0.5
    dispersivity_m = 0.5
    
    # Derived environmental parameters
    # The term 'Cgw' seems to represent the bulk density of the porous medium.
    # Convert density to kg/L for unit consistency with Kd (L/kg).
    Cgw_density_kg_L = coal_seam_density_kg_m3 / 1000.0
    # Dispersion Coefficient D = advection velocity * dispersivity factor
    D_dispersion_m2_d = v_velocity_m_d * dispersivity_m
    
    # Chemical-specific parameters in a dictionary
    chemicals = {
        'Atrazine': {
            'product_fraction': 0.01,
            'conc_in_product_ug_L': 40.0,
            'log_koc': 2.20,
            'half_life_d': 90.0,
            'ec50_ug_L': 100.0,
        },
        'PFOS': {
            'product_fraction': 0.001,
            'conc_in_product_ug_L': 300.0,
            'log_koc': 3.65,
            'half_life_d': 14965.0,
            'ec50_ug_L': 480.0,
        },
        'Endosulfan': {
            'product_fraction': 0.005,
            'conc_in_product_ug_L': 20.0,
            'log_koc': 4.3,
            'half_life_d': 60.0,
            'ec50_ug_L': 560.0,
        }
    }
    
    # --- Step 2: Calculate derived parameters for each chemical ---
    for name, params in chemicals.items():
        # MTotal: Total mass of chemical injected (μg)
        params['MTotal_ug'] = total_water_volume_L * params['product_fraction'] * params['conc_in_product_ug_L']
        
        # k: First-order decay rate constant (1/d)
        params['k_decay_rate'] = math.log(2) / params['half_life_d']
        
        # Kd: Partitioning coefficient (L/kg) = Koc * foc
        koc = 10**params['log_koc']
        params['Kd_L_kg'] = koc * coal_foc
        
        # Retardation term from the user's formula
        params['retardation_term'] = 1.0 / (1.0 + params['Kd_L_kg'] * Cgw_density_kg_L)

    def calculate_concentration(t, params):
        """Calculates concentration C(x,t) using the provided formula."""
        if t <= 0:
            return 0.0
        
        MTotal = params['MTotal_ug']
        k = params['k_decay_rate']
        retard_term = params['retardation_term']
        
        # Denominator parts
        sqrt_term = math.sqrt(4 * math.pi * D_dispersion_m2_d * t)
        denom_exp_term = 4 * D_dispersion_m2_d * t

        # Exponent parts
        num_exp_term = -((x_distance_m - v_velocity_m_d * t)**2)
        
        # Handle potential math domain errors for exp
        try:
            exp_adv_disp = math.exp(num_exp_term / denom_exp_term)
        except OverflowError:
            exp_adv_disp = 0 # exp(-infinity) is 0

        exp_decay = math.exp(-k * t)
        
        # Full formula
        concentration = (MTotal / sqrt_term) * exp_adv_disp * exp_decay * retard_term
        return concentration

    # --- Step 3 & 4: Find the peak concentration for each chemical ---
    peak_info = {}
    time_range = range(1, 1001) # Search for peak time over a reasonable range of days

    for name, params in chemicals.items():
        max_conc = 0
        peak_time = 0
        for t_day in time_range:
            conc = calculate_concentration(t_day, params)
            if conc > max_conc:
                max_conc = conc
                peak_time = t_day
        peak_info[name] = {'peak_conc_ug_L': max_conc, 'peak_time_d': peak_time}

    # Find the overall highest concentration
    highest_conc = 0
    highest_chem = None
    time_at_highest_conc = 0

    for name, info in peak_info.items():
        if info['peak_conc_ug_L'] > highest_conc:
            highest_conc = info['peak_conc_ug_L']
            highest_chem = name
            time_at_highest_conc = info['peak_time_d']

    print(f"The highest concentration of an individual chemical reaching the spring is {highest_conc:.2f} μg/L.")
    print(f"This is from {highest_chem}, occurring at day {time_at_highest_conc}.")
    print("-" * 30)

    # --- Step 5: Assess mixture toxicity at the time of the peak concentration ---
    print(f"Calculating mixture effect at the time of the peak: {time_at_highest_conc} days.\n")

    # Calculate concentration of all chemicals at this specific time
    concentrations_at_peak_time = {}
    for name, params in chemicals.items():
        concentrations_at_peak_time[name] = calculate_concentration(time_at_highest_conc, params)

    # Calculate Hazard Index (HI)
    HI = 0
    hi_equation_str = "Hazard Index = "
    
    for name, conc in concentrations_at_peak_time.items():
        ec50 = chemicals[name]['ec50_ug_L']
        tu = conc / ec50 # Toxicity Unit
        HI += tu
        print(f"Concentration of {name} at day {time_at_highest_conc}: {conc:.4f} μg/L")
        hi_equation_str += f"({conc:.4f} / {ec50}) + "

    # --- Step 6: Determine the effect and print the final equation ---
    print("\nThe Hazard Index (HI) is calculated as the sum of concentrations divided by their EC50 values:")
    # Print the final equation with values
    print(f"{hi_equation_str.strip(' + ')}")
    print(f"Hazard Index = {HI:.4f}")

    if HI > 1.1:
        effect = "Synergistic (or greater than additive)"
    elif HI < 0.9:
        effect = "Antagonistic (or less than additive)"
    else:
        effect = "Additive"
        
    print(f"\nBased on the Hazard Index of {HI:.4f}, the mixture effect is considered to be {effect}.")
    
    # Return the final answer in the required format
    return f"{highest_conc:.2f} μg/L, {effect}"

# Run the simulation and print the final answer in the specified format
final_answer = solve_environmental_problem()
print(f"\n<<<C({final_answer})>>>")