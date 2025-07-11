import math

def solve_contaminant_transport():
    """
    Calculates the highest chemical concentration reaching a spring and the mixture effect
    based on a simplified advection-dispersion model.
    """
    # 1. Define constants and parameters from the problem
    params = {
        'env': {
            'injected_water_L': 1500000,
            'coal_density_kg_m3': 1346,
            'f_oc': 0.50,
            'porosity_n': 0.30,
            'fracture_height_m': 10,
            'fracture_width_m': 0.005,
            'distance_to_spring_m': 100,
            'advection_velocity_m_d': 0.5,
            'dispersivity_m': 0.5,
        },
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
            'log_koc': 4.3,
            'half_life_d': 60,
            'ec50_ug_L': 560,
        }
    }

    # Environmental parameters
    env = params['env']
    V_water_m3 = env['injected_water_L'] / 1000.0
    x = env['distance_to_spring_m']
    v = env['advection_velocity_m_d']
    n = env['porosity_n']
    rho_b = env['coal_density_kg_m3']
    f_oc = env['f_oc']
    A = env['fracture_height_m'] * env['fracture_width_m']
    D = env['dispersivity_m'] * v
    
    # Simplified model: peak travels with groundwater velocity
    t_peak = x / v

    results = {}
    
    # 2. Loop through each chemical to calculate parameters and concentrations
    for name, chem in params.items():
        if name == 'env':
            continue

        # Mass (M) in kg
        # ug/L is equivalent to mg/m3. Convert to kg/m3 by dividing by 1000.
        conc_in_product_kg_m3 = chem['conc_in_product_ug_L'] / 1000.0
        V_product_m3 = V_water_m3 * chem['product_percent']
        M = V_product_m3 * conc_in_product_kg_m3

        # Decay constant (k) in d^-1
        k = math.log(2) / chem['half_life_d']

        # Retardation factor (R)
        Koc = 10**chem['log_koc']  # L/kg
        Kd_L_kg = Koc * f_oc
        Kd_m3_kg = Kd_L_kg / 1000.0
        R = 1 + (rho_b / n) * Kd_m3_kg

        # 3. Calculate peak concentration using the simplified model
        # C_peak = (M / (A*n)) * (1 / sqrt(4*pi*D*t)) * exp(-k*t) * (1/R)
        # This is equivalent to C_peak = (M / (A*n*R)) * (1/sqrt(...)) * exp(...)
        
        source_term = M / (A * n * R)
        dispersion_term = 1 / math.sqrt(4 * math.pi * D * t_peak)
        decay_term = math.exp(-k * t_peak)
        
        C_peak_kg_m3 = source_term * dispersion_term * decay_term
        C_peak_ug_L = C_peak_kg_m3 * 1e6 # 1 kg/m3 = 1e6 ug/L

        results[name] = {
            'M_kg': M,
            'k_d-1': k,
            'R': R,
            'C_peak_ug_L': C_peak_ug_L,
            'EC50': chem['ec50_ug_L']
        }

    # 4. Find the highest concentration
    highest_conc = 0
    highest_conc_chem = None
    for name, data in results.items():
        if data['C_peak_ug_L'] > highest_conc:
            highest_conc = data['C_peak_ug_L']
            highest_conc_chem = name

    print("--- Part 1: Highest Individual Chemical Concentration ---")
    print(f"Based on a simplified transport model, the peak of the contaminant plume reaches the spring at t = {t_peak:.1f} days.")
    print(f"The chemical with the highest concentration at this time is: {highest_conc_chem}\n")
    
    # Print the detailed calculation for the winning chemical
    winner_data = results[highest_conc_chem]
    print("Calculation for highest concentration (Atrazine):")
    print("C_peak = (M / (A * n * R)) * (1 / sqrt(4 * pi * D * t_peak)) * exp(-k * t_peak)")
    print(f"C_peak = ({winner_data['M_kg']:.4f} kg / ({A:.3f} m^2 * {n:.2f} * {winner_data['R']:.2f})) * "
          f"(1 / sqrt(4 * pi * {D:.2f} m^2/d * {t_peak:.1f} d)) * "
          f"exp(-{winner_data['k_d-1']:.4f} d^-1 * {t_peak:.1f} d)")
    print(f"Highest Concentration = {highest_conc:.3f} μg/L\n")

    # 5. Assess mixture toxicity
    print("--- Part 2: Mixture Effect Analysis ---")
    print(f"The concentrations of all chemicals at t = {t_peak:.1f} days are:")
    for name, data in results.items():
        print(f"  - C_{name}: {data['C_peak_ug_L']:.3f} μg/L")

    print("\nThe mixture effect is assessed using the Toxic Unit (TU) model, which assumes additive effects.")
    print("TU_mix = sum(Concentration_i / EC50_i)")
    
    tu_mix = 0
    tu_calcs = []
    for name, data in results.items():
        tu = data['C_peak_ug_L'] / data['EC50']
        tu_mix += tu
        tu_calcs.append(f"({data['C_peak_ug_L']:.3f} / {data['EC50']})")
    
    print(f"TU_mix = {' + '.join(tu_calcs)}")
    
    tu_values = [res['C_peak_ug_L'] / res['EC50'] for res in results.values()]
    print(f"TU_mix = {tu_values[0]:.5f} + {tu_values[1]:.5f} + {tu_values[2]:.5f}")
    print(f"TU_mix = {tu_mix:.5f}\n")

    print("Conclusion:")
    print("The Mixture Toxic Unit is significantly less than 1, indicating a low toxic risk to the algal community from the mixture at these concentrations.")
    print("The toxic unit model assumes an ADDITIVE effect, where the combined toxicity is the sum of the individual toxicities scaled by their potency.")

    # Final answer in the required format
    final_answer = f"The highest concentration is {highest_conc:.3f} μg/L (Atrazine), and the assumed mixture effect is Additive."
    return final_answer

# Run the analysis and print the final answer
final_answer_string = solve_contaminant_transport()
print(f"\n<<<{final_answer_string}>>>")
