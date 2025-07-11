import math

def solve_environmental_risk():
    """
    This function solves the contaminant transport and risk problem as described.
    It calculates the highest concentration of three chemicals reaching a spring
    and determines the nature of the mixture's effect.
    """

    # Step 1: Define constants from the problem description
    print("1. Defining constants and environmental parameters...")
    # Environmental parameters
    x_distance = 100.0  # m, distance to the spring
    v_velocity = 0.5  # m/d, advection velocity
    alpha_dispersivity = 0.5  # m, dispersivity factor
    D_dispersion = alpha_dispersivity * v_velocity  # m^2/d, dispersion coefficient
    
    # Coal seam properties
    rho_b_density_m3 = 1346.0  # kg/m3, coal seam density
    rho_b_density_L = rho_b_density_m3 / 1000.0 # kg/L, for unit consistency with Kd
    foc_organic_carbon = 0.50  # fraction of organic carbon
    
    # Fracture (source) properties
    fracture_length = 80.0  # m
    fracture_height = 10.0  # m
    source_area = fracture_length * fracture_height  # m^2

    # Injection properties
    total_water_volume_L = 1500000.0 # L

    # Step 2: Define chemical-specific properties
    chemicals = {
        'Atrazine': {
            'product_percentage': 0.01,
            'product_conc_ug_L': 40.0,
            'log_Koc': 2.20,
            'half_life_days': 90.0,
            'EC50_ug_L': 100.0,
        },
        'PFOS': {
            'product_percentage': 0.001,
            'product_conc_ug_L': 300.0,
            'log_Koc': 3.65,
            'half_life_days': 14965.0,
            'EC50_ug_L': 480.0,
        },
        'Endosulfan': {
            'product_percentage': 0.005,
            'product_conc_ug_L': 20.0,
            'log_Koc': 4.3,
            'half_life_days': 60.0,
            'EC50_ug_L': 560.0,
        }
    }
    print("...Done.\n")

    # Step 3 & 4: Perform calculations for each chemical to find peak concentration
    print("2. Calculating peak concentration for each chemical...")
    results = {}
    
    for name, params in chemicals.items():
        # a. Calculate total mass and mass per unit area (M_area)
        volume_of_product = total_water_volume_L * params['product_percentage']
        total_mass_ug = volume_of_product * params['product_conc_ug_L']
        total_mass_g = total_mass_ug / 1e6
        M_area = total_mass_g / source_area  # g/m^2
        params['M_area'] = M_area

        # b. Calculate decay rate constant (k)
        k = math.log(2) / params['half_life_days']  # 1/d
        params['k'] = k
        
        # c. Calculate the sorption factor from the term '1/(1+KdCgw)'
        Koc = 10**params['log_Koc']
        Kd = Koc * foc_organic_carbon # L/kg
        sorption_factor = 1.0 / (1.0 + Kd * rho_b_density_L)
        params['sorption_factor'] = sorption_factor
        
        # d. Calculate t_peak (time of peak concentration) by solving the quadratic equation:
        # (v^2 + 4Dk)t^2 + 2Dt - x^2 = 0
        a_quad = v_velocity**2 + 4 * D_dispersion * k
        b_quad = 2 * D_dispersion
        c_quad = -x_distance**2
        t_peak = (-b_quad + math.sqrt(b_quad**2 - 4 * a_quad * c_quad)) / (2 * a_quad)
        params['t_peak'] = t_peak

        # e. Calculate C_peak (peak concentration) using the provided formula
        sqrt_term = math.sqrt(4 * math.pi * D_dispersion * t_peak)
        exp_advection_num = -(x_distance - v_velocity * t_peak)**2
        exp_advection_den = 4 * D_dispersion * t_peak
        exp_advection_term = math.exp(exp_advection_num / exp_advection_den)
        exp_decay_term = math.exp(-k * t_peak)
        
        C_peak_g_m3 = (M_area / sqrt_term) * exp_advection_term * exp_decay_term * sorption_factor
        C_peak_ug_L = C_peak_g_m3 * 1000.0 # Conversion from g/m^3 to ug/L
        
        results[name] = {
            'peak_conc_ug_L': C_peak_ug_L,
            'peak_time_days': t_peak,
            'params': params
        }
    print("...Done.\n")
    
    # Step 5: Find the highest concentration and display the full equation
    highest_conc = -1.0
    highest_conc_chem = None
    for name, res in results.items():
        if res['peak_conc_ug_L'] > highest_conc:
            highest_conc = res['peak_conc_ug_L']
            highest_conc_chem = name
            
    print(f"3. Identifying the highest concentration and showing the calculation...")
    print(f"The highest concentration is from {highest_conc_chem}.\n")
    
    chem_details = results[highest_conc_chem]
    p = chem_details['params']
    t_p = chem_details['peak_time_days']

    print(f"The full calculation for {highest_conc_chem} is:")
    print("C(x,t) = (M_area / sqrt(4*pi*D*t)) * exp(-(x-vt)^2/(4Dt)) * exp(-kt) * SorptionFactor\n")
    print(f"C({x_distance:.0f}, {t_p:.2f}) = ({p['M_area']:.4e} g/m^2 / sqrt(4 * pi * {D_dispersion:.2f} m^2/d * {t_p:.2f} d)) * \\")
    print(f"    exp(-({x_distance:.0f} m - {v_velocity:.1f} m/d * {t_p:.2f} d)^2 / (4 * {D_dispersion:.2f} m^2/d * {t_p:.2f} d)) * \\")
    print(f"    exp(-{p['k']:.4e} 1/d * {t_p:.2f} d) * {p['sorption_factor']:.4e}\n")
    
    # Step 6: State the final answers
    mixture_effect_type = "additive"
    
    print("-------------------- FINAL ANSWER --------------------")
    print(f"The highest concentration of an individual chemical reaching the spring is {highest_conc:.4e} μg/L.")
    print(f"The effect of the mixture is assumed to be {mixture_effect_type}, based on standard risk assessment principles for chemical mixtures without specified interaction data.")
    print("----------------------------------------------------")
    
    final_output_val = f"{highest_conc:.4e} μg/L, {mixture_effect_type}"
    print(f"\n<<<The highest concentration is {highest_conc:.4e} ug/L and the mixture effect is additive.>>>")


# Execute the main function
solve_environmental_risk()