import math

def solve_transport():
    # --- Step 1: Define Constants and Initial Conditions ---
    # Environmental parameters
    x = 100.0  # distance to spring (m)
    v = 0.5  # advection velocity (m/d)
    dispersivity = 0.5  # dispersivity factor (m)
    D = v * dispersivity  # Dispersion coefficient (m^2/d)
    
    # Coal seam properties
    coal_density = 1346.0  # kg/m^3
    foc = 0.50  # fraction of organic carbon
    porosity = 0.30  # water content
    
    # Fracture and injection properties
    fracture_length = 80.0  # m
    fracture_height = 10.0  # m
    source_area = fracture_length * fracture_height # m^2
    total_water_injected = 1500000.0 # L

    # --- Step 2: Define Chemical Properties ---
    chemicals = {
        'Atrazine': {
            'product_percent': 0.01,
            'conc_in_product': 40.0,  # ug/L
            'log_koc': 2.20,
            'half_life': 90.0,  # days
            'ec50': 100.0  # ug/L
        },
        'PFOS': {
            'product_percent': 0.001,
            'conc_in_product': 300.0, # ug/L
            'log_koc': 3.65,
            'half_life': 14965.0, # days
            'ec50': 480.0  # ug/L
        },
        'Endosulfan': {
            'product_percent': 0.005,
            'conc_in_product': 20.0,  # ug/L
            'log_koc': 4.3,
            'half_life': 60.0,  # days
            'ec50': 560.0  # ug/L
        }
    }
    
    results = {}
    total_hazard_index = 0

    # --- Step 3: Loop through each chemical and perform calculations ---
    for name, props in chemicals.items():
        # Calculate total mass injected (M) in ug
        mass_total = total_water_injected * props['product_percent'] * props['conc_in_product']
        
        # Calculate MTotal for the formula (mass per unit area, M/A) in ug/m^2
        m_per_area = mass_total / source_area
        
        # Calculate decay rate constant (k) in 1/d
        k = math.log(2) / props['half_life']
        
        # Calculate partition coefficient (Kd) and Retardation Factor (R)
        koc = 10**props['log_koc']
        kd = koc * foc  # L/kg
        # R = 1 + (bulk_density / porosity) * Kd. We convert Kd from L/kg to m3/kg by dividing by 1000
        R = 1 + (coal_density / porosity) * (kd / 1000.0)
        
        # Find time of peak concentration (t_peak) by solving quadratic equation At^2 + Bt + C = 0
        quad_a = v**2 + 4 * D * k
        quad_b = 2 * D
        quad_c = -x**2
        
        # Using quadratic formula: t = (-b + sqrt(b^2 - 4ac)) / 2a
        t_peak = (-quad_b + math.sqrt(quad_b**2 - 4 * quad_a * quad_c)) / (2 * quad_a)
        
        # --- Step 4: Calculate the Maximum Concentration (C_max) ---
        # Calculate individual terms of the equation for clarity
        exp_term_1 = math.exp(-((x - v * t_peak)**2) / (4 * D * t_peak))
        exp_term_2 = math.exp(-k * t_peak)
        sqrt_term = math.sqrt(4 * math.pi * D * t_peak)
        
        # Calculate C_max in ug/m^3 based on the interpreted formula
        c_max_ug_m3 = (m_per_area / sqrt_term) * exp_term_1 * exp_term_2 * (1/R)
        
        # Convert C_max to ug/L (1 m^3 = 1000 L)
        c_max_ug_L = c_max_ug_m3 / 1000.0
        
        # --- Step 5: Calculate Hazard Quotient and store results ---
        hazard_quotient = c_max_ug_L / props['ec50']
        total_hazard_index += hazard_quotient
        
        results[name] = {
            'c_max_ug_L': c_max_ug_L,
            'hazard_quotient': hazard_quotient,
            't_peak': t_peak,
            'm_per_area': m_per_area,
            'k': k,
            'R': R
        }

    # --- Step 6: Identify highest concentration and present results ---
    highest_conc_chem = max(results, key=lambda k: results[k]['c_max_ug_L'])
    highest_conc_val = results[highest_conc_chem]['c_max_ug_L']

    print(f"Calculation for the chemical with the highest concentration: {highest_conc_chem}\n")

    # Retrieve values for printing the equation
    chem_props = results[highest_conc_chem]
    m_val = chem_props['m_per_area']
    d_val = D
    t_val = chem_props['t_peak']
    x_val = x
    v_val = v
    k_val = chem_props['k']
    r_val = chem_props['R']

    print("The final concentration is calculated using the interpreted formula:")
    print("C(x,t) = (M/A) * (1/sqrt(4*pi*D*t)) * exp(-(x-v*t)^2/(4*D*t)) * exp(-k*t) * (1/R)\n")

    print("Plugging in the numbers for the peak concentration of Atrazine:")
    print(f"C_max [ug/L] = (1/1000) * ( {m_val:.2f} / sqrt(4 * pi * {d_val:.2f} * {t_val:.2f}) * exp(-({x_val:.2f} - {v_val:.2f}*{t_val:.2f})^2 / (4 * {d_val:.2f} * {t_val:.2f})) * exp(-{k_val:.5f} * {t_val:.2f}) * (1 / {r_val:.2f}) )")
    print(f"\nThe highest concentration reached at the spring is {highest_conc_val:.3e} ug/L from {highest_conc_chem}.")

    print("\n--- Mixture Effect Analysis ---")
    print(f"The Hazard Index (HI) for the mixture is the sum of individual Hazard Quotients (HQ = C_max/EC50).")
    print(f"HI = {results['Atrazine']['hazard_quotient']:.2e} (Atrazine) + {results['PFOS']['hazard_quotient']:.2e} (PFOS) + {results['Endosulfan']['hazard_quotient']:.2e} (Endosulfan)")
    print(f"Total HI = {total_hazard_index:.2e}")

    if total_hazard_index < 1:
        print("Since the Hazard Index is much less than 1, the risk from the mixture is considered negligible.")
    else:
        print("Since the Hazard Index is greater than or equal to 1, there is a potential risk from the mixture.")
        
    print("\nBased on the Hazard Index calculation, the mixture effect is assumed to be additive, as there is no data to suggest synergistic or antagonistic effects.")
    
    # Final answer for the system
    final_answer_str = f"The highest concentration is {highest_conc_val:.3e} μg/L from {highest_conc_chem}. The mixture effect is assumed to be additive."
    print(f"\n<<<{final_answer_str}>>>")


solve_transport()