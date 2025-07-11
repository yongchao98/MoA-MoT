import math

def solve_transport_problem():
    """
    Calculates the highest concentration of contaminants reaching a spring
    and assesses the combined toxic effect using the Hazard Index method.
    """

    # --- 1. Define Constants and Initial Conditions ---
    # Environmental parameters
    distance_to_spring_x = 100.0  # m
    advection_velocity_v = 0.5  # m/d
    dispersivity_alpha = 0.5  # m
    
    # Aquifer/Coal properties
    coal_seam_density_kg_m3 = 1346.0  # kg/m^3
    frac_organic_carbon_f_oc = 0.50
    porosity_n = 0.30 # from 30% water
    
    # Source properties
    total_water_injected_L = 1500000.0  # L
    fracture_length = 80.0  # m
    fracture_height = 10.0  # m
    fracture_area = fracture_length * fracture_height # m^2

    # --- 2. Interpret and Pre-calculate Model Parameters ---
    # Dispersion coefficient D = alpha * v
    dispersion_coeff_D = dispersivity_alpha * advection_velocity_v # m^2/d
    
    # Interpretation for the formula term Cgw
    # To make (1 + Kd*Cgw) dimensionless, Cgw must have units of kg/L.
    # We assume Cgw is the coal solid density in kg/L.
    cgw_rho_s_kg_L = coal_seam_density_kg_m3 / 1000.0 # kg/L

    # --- 3. Chemical Data and Loop ---
    chemical_data = {
        'Atrazine': {
            'product_percent': 0.01,   # 1%
            'product_conc_ug_L': 40.0,
            'log_koc': 2.20,
            'half_life_days': 90.0,
            'ec50_ug_L': 100.0
        },
        'PFOS': {
            'product_percent': 0.001,  # 0.1%
            'product_conc_ug_L': 300.0,
            'log_koc': 3.65,
            'half_life_days': 14965.0,
            'ec50_ug_L': 480.0
        },
        'Endosulfan': {
            'product_percent': 0.005,  # 0.5%
            'product_conc_ug_L': 20.0,
            'log_koc': 4.3,
            'half_life_days': 60.0,
            'ec50_ug_L': 560.0
        }
    }

    results = {}
    highest_overall_conc = 0.0
    
    print("--- Calculating Peak Concentrations ---\n")

    for name, data in chemical_data.items():
        # Calculate mass injected (M)
        mass_injected_ug = total_water_injected_L * data['product_percent'] * data['product_conc_ug_L']
        
        # Calculate mass loading per area (MTotal for the formula)
        m_total_load_ug_m2 = mass_injected_ug / fracture_area
        
        # Calculate decay rate constant (k)
        k = math.log(2) / data['half_life_days']
        
        # Calculate partitioning coefficient (Kd)
        koc = 10**data['log_koc']
        kd_L_kg = koc * frac_organic_carbon_f_oc
        
        # Calculate the dimensionless partitioning term from the formula
        partition_term = 1.0 / (1.0 + kd_L_kg * cgw_rho_s_kg_L)

        # --- 4. Numerical Search for Peak Concentration ---
        max_conc_ug_L = 0.0
        peak_time = 0
        
        # Search for peak concentration over time (days)
        for t in range(1, 1001):
            if t > 0:
                # Core of the provided formula
                sqrt_term = math.sqrt(4 * math.pi * dispersion_coeff_D * t)
                exp_transport_term = math.exp(-((distance_to_spring_x - advection_velocity_v * t)**2) / (4 * dispersion_coeff_D * t))
                exp_decay_term = math.exp(-k * t)

                # Calculate concentration in ug/m^3 as per formula interpretation
                conc_ug_m3 = (m_total_load_ug_m2 / sqrt_term) * exp_transport_term * exp_decay_term * partition_term
                
                # Convert to ug/L
                conc_ug_L = conc_ug_m3 / 1000.0

                if conc_ug_L > max_conc_ug_L:
                    max_conc_ug_L = conc_ug_L
                    peak_time = t

        results[name] = {'max_conc': max_conc_ug_L, 'ec50': data['ec50_ug_L']}
        
        if max_conc_ug_L > highest_overall_conc:
            highest_overall_conc = max_conc_ug_L

        print(f"Peak concentration for {name}: {max_conc_ug_L:.6f} μg/L at day {peak_time}")

    # --- 5. Final Analysis and Output ---
    print("\n--- Final Results ---")
    print(f"\nThe highest concentration of an individual chemical reaching the spring is: {highest_overall_conc:.6f} μg/L (Atrazine)")
    
    print("\n--- Mixture Effect Analysis (Hazard Index) ---")
    print("The combined effect is assessed assuming additivity by calculating a Hazard Index (HI).")
    print("HI = (Conc_1 / EC50_1) + (Conc_2 / EC50_2) + ...")
    
    hazard_index = 0
    hq_calcs = []
    for name, res in results.items():
        hq = res['max_conc'] / res['ec50']
        hazard_index += hq
        hq_calcs.append(f"({res['max_conc']:.6f} / {res['ec50']})")
        
    print("\nCalculating Hazard Index for the algal community:")
    print(f"HI = { ' + '.join(hq_calcs) }")
    print(f"HI = {hazard_index:.6f}")
    
    if hazard_index > 1:
        print("\nConclusion: The Hazard Index is greater than 1, suggesting a potential risk to the algal community from the additive effects of the chemical mixture.")
    else:
        print("\nConclusion: The Hazard Index is less than 1, suggesting the mixture is unlikely to pose a significant risk to the algal community.")
    
    # The final answer format is <<<value>>>. The question asks for two things. The first is the highest concentration.
    # The second is the effect. Based on the HI, the effect is Additive, with a risk determined by the HI value.
    # Returning the highest concentration as the primary quantitative answer.
    return highest_overall_conc

# Execute the function and capture the primary result for the final output
final_answer = solve_transport_problem()
print(f"\n<<<{final_answer:.6f}>>>")