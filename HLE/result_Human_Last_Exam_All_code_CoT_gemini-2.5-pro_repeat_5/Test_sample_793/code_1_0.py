import math

def solve_concentration_and_effect():
    """
    Calculates the peak concentration of contaminants reaching a spring and assesses the mixture effect.
    """
    # --- Step 1: Deconstruct and Parameterize ---
    # General Parameters
    total_water_L = 1_500_000
    dist_x = 100.0  # m
    velo_v = 0.5  # m/d
    dispersivity_alpha = 0.5  # m
    fracture_L = 80.0 # m
    fracture_H = 10.0 # m
    
    # Aquifer/Coal Properties
    coal_density_kg_m3 = 1346.0
    foc = 0.50  # fraction of organic carbon
    water_content_theta = 0.30  # porosity

    # Chemical Properties in a list of dictionaries
    chemicals = [
        {
            'name': 'Atrazine',
            'product_percent': 0.01,
            'conc_ug_L': 40.0,
            'log_koc': 2.20,
            'half_life_d': 90.0,
            'ec50_ug_L': 100.0,
        },
        {
            'name': 'PFOS',
            'product_percent': 0.001,
            'conc_ug_L': 300.0,
            'log_koc': 3.65,
            'half_life_d': 14965.0,
            'ec50_ug_L': 480.0,
        },
        {
            'name': 'Endosulfan',
            'product_percent': 0.005,
            'conc_ug_L': 20.0,
            'log_koc': 4.3,
            'half_life_d': 60.0,
            'ec50_ug_L': 560.0,
        },
    ]

    results = []
    total_hq = 0.0

    print("Calculating contaminant transport and effects...\n")

    for chem in chemicals:
        # --- Step 2: Calculate Intermediate Values ---

        # Total mass injected (ug)
        m_total_ug = total_water_L * chem['product_percent'] * chem['conc_ug_L']

        # Source Area (m^2) and Mass per Area (ug/m^2)
        source_area_m2 = fracture_L * fracture_H
        m_area_ug_m2 = m_total_ug / source_area_m2

        # Dispersion coefficient (m^2/d)
        D_m2_d = dispersivity_alpha * velo_v

        # Decay rate constant (1/d)
        k_per_d = math.log(2) / chem['half_life_d']

        # Partitioning coefficient Kd (L/kg) and Retardation Factor term
        koc = 10**chem['log_koc']
        kd_L_kg = koc * foc
        
        # Cgw term interpreted as (bulk density in kg/L) / porosity
        cgw_term_kg_L = (coal_density_kg_m3 / 1000.0) / water_content_theta
        partition_factor = 1.0 / (1.0 + kd_L_kg * cgw_term_kg_L)

        # --- Step 3: Find Time of Peak Concentration (t_peak) ---
        # Solves quadratic equation: (v^2 + 4Dk)t^2 + 2Dt - x^2 = 0
        a = velo_v**2 + 4 * D_m2_d * k_per_d
        b = 2 * D_m2_d
        c = -dist_x**2
        
        # Using the quadratic formula, taking the positive root for time
        t_peak_d = (-b + math.sqrt(b**2 - 4 * a * c)) / (2 * a)

        # --- Step 4: Calculate Peak Concentration (C_max) ---
        t = t_peak_d
        
        # Main equation terms
        term1_sqrt = math.sqrt(4 * math.pi * D_m2_d * t)
        
        exp_num = -(dist_x - velo_v * t)**2
        exp_den = 4 * D_m2_d * t
        term2_exp_transport = math.exp(exp_num / exp_den)
        
        term3_exp_decay = math.exp(-k_per_d * t)

        # Calculate concentration in ug/m^3
        c_max_ug_m3 = (m_area_ug_m2 / term1_sqrt) * term2_exp_transport * term3_exp_decay * partition_factor
        
        # Convert to ug/L
        c_max_ug_L = c_max_ug_m3 / 1000.0

        # --- Step 5: Assess Mixture Effect ---
        hq = c_max_ug_L / chem['ec50_ug_L']
        total_hq += hq
        
        # Store results
        results.append({
            'name': chem['name'],
            'C_max_ug_L': c_max_ug_L,
            'HQ': hq,
            't_peak_d': t_peak_d,
            'm_area_ug_m2': m_area_ug_m2,
            'D_m2_d': D_m2_d,
            'k_per_d': k_per_d,
            'kd_L_kg': kd_L_kg,
            'cgw_term_kg_L': cgw_term_kg_L,
            'partition_factor': partition_factor
        })

    # --- Step 6: Output Results ---
    highest_conc_chem = max(results, key=lambda x: x['C_max_ug_L'])

    print("--- Individual Chemical Results ---\n")
    for res in results:
        print(f"Chemical: {res['name']}")
        print(f"  Peak Concentration (C_max): {res['C_max_ug_L']:.4e} ug/L")
        print(f"  Time to Peak: {res['t_peak_d']:.2f} days")
        print(f"  Hazard Quotient (HQ): {res['HQ']:.4e}\n")

    print("\n--- Summary and Conclusion ---\n")
    print(f"The chemical with the highest concentration reaching the spring is {highest_conc_chem['name']}.")
    
    print("\nThe full equation for this chemical is:")
    print("C(x,t) = [ (M_area / sqrt(4*pi*D*t)) * exp(-(x-v*t)^2/(4*D*t)) * exp(-k*t) * (1/(1+Kd*Cgw)) ] / 1000")
    
    # Print the equation with the final numbers
    res = highest_conc_chem
    print("C(100, {t:.2f}) = [ ({m_area:.2f} / sqrt(4*pi*{D:.2f}*{t:.2f})) * exp(-({x:.2f}-{v:.2f}*{t:.2f})^2 / (4*{D:.2f}*{t:.2f})) * exp(-{k:.4f}*{t:.2f}) * {pf:.6f} ] / 1000 = {c_max:.4e} ug/L".format(
        t=res['t_peak_d'],
        m_area=res['m_area_ug_m2'],
        D=res['D_m2_d'],
        x=dist_x,
        v=velo_v,
        k=res['k_per_d'],
        pf=res['partition_factor'],
        c_max=res['C_max_ug_L']
    ))
    
    print(f"\nMixture Effect Assessment:")
    print(f"The total Hazard Index (HI) for the mixture is {total_hq:.4e}.")
    print("Since the Hazard Index (HI) is significantly less than 1, the risk to the algal community from the mixture is considered negligible.")
    print("The Hazard Index method assumes an ADDITIVE effect, where the combined toxicity is the sum of the individual chemical toxicities scaled by their concentration.")

    # Final answer format
    final_answer_str = f"The highest concentration is {highest_conc_chem['C_max_ug_L']:.4e} ug/L from {highest_conc_chem['name']}. The effect of the mixture is assessed using an additive model, which results in a Hazard Index of {total_hq:.4e}, indicating a negligible effect."
    print(f"\n<<<{final_answer_str}>>>")


solve_concentration_and_effect()