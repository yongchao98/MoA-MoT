import math

def solve_transport_problem():
    """
    Calculates the highest concentration of chemicals reaching a spring
    and assesses the mixture effect based on a given advection-diffusion equation.
    """
    # --- 1. Define Initial Parameters ---
    # Physical and transport parameters
    injected_water_L = 1500000
    coal_density_kg_m3 = 1346
    frac_organic_carbon = 0.50
    porosity = 0.30
    fracture_length_m = 80
    fracture_height_m = 10
    distance_to_spring_m = 100
    advection_velocity_m_d = 0.5
    dispersivity_m = 0.5

    # Chemical properties
    chemicals_data = {
        'Atrazine': {
            'product_percent': 0.01,
            'conc_in_product_ug_L': 40,
            'log_koc': 2.20,
            'half_life_d': 90,
            'ec50_ug_L': 100
        },
        'PFOS': {
            'product_percent': 0.001,
            'conc_in_product_ug_L': 300,
            'log_koc': 3.65,
            'half_life_d': 14965,
            'ec50_ug_L': 480
        },
        'Endosulfan': {
            'product_percent': 0.005,
            'conc_in_product_ug_L': 20,
            'log_koc': 4.3,
            'half_life_d': 60,
            'ec50_ug_L': 560
        }
    }
    
    # --- Derived parameters ---
    dispersion_coeff_D = advection_velocity_m_d * dispersivity_m # m^2/d
    fracture_area_m2 = fracture_length_m * fracture_height_m
    x = distance_to_spring_m
    v = advection_velocity_m_d
    D = dispersion_coeff_D
    rho_b = coal_density_kg_m3
    theta = porosity

    print("--- Problem Setup ---")
    print(f"Distance to spring (x): {x} m")
    print(f"Groundwater velocity (v): {v} m/d")
    print(f"Dispersion coefficient (D): {D} m^2/d")
    print(f"Fracture face area for mass loading: {fracture_area_m2} m^2")
    print("-" * 25)

    results = []
    
    # --- 2. Calculations for each chemical ---
    for name, props in chemicals_data.items():
        print(f"\nCalculating for {name}...")
        
        # Total mass in grams
        mass_ug = injected_water_L * props['product_percent'] * props['conc_in_product_ug_L']
        mass_g = mass_ug / 1e6
        
        # Mass loading MTotal (M/A) in g/m^2
        MTotal = mass_g / fracture_area_m2
        print(f"  - Total Injected Mass: {mass_g:.4f} g")
        print(f"  - Mass Loading (MTotal): {MTotal:.3e} g/m^2")
        
        # Decay constant k in 1/d
        k = math.log(2) / props['half_life_d']
        print(f"  - Decay Constant (k): {k:.4f} 1/d")
        
        # Partitioning coefficient Kd
        koc = 10**props['log_koc']  # L/kg
        koc_m3_kg = koc * 0.001  # m^3/kg
        Kd = frac_organic_carbon * koc_m3_kg # m^3/kg
        print(f"  - Partition Coefficient (Kd): {Kd:.4f} m^3/kg")
        
        # Retardation-like factor F from 1/(1 + Kd*Cgw) = 1/(1 + Kd*(rho_b/theta))
        F = 1 / (1 + Kd * (rho_b / theta))
        print(f"  - Retardation-like Factor (F): {F:.3e}")
        
        # Time of peak concentration t_peak (quadratic solution: a*t^2 + b*t + c = 0)
        a = v**2 + 4 * D * k
        b = 2 * D
        c = -x**2
        t_peak = (-b + math.sqrt(b**2 - 4 * a * c)) / (2 * a)
        print(f"  - Time of Peak Concentration (t_peak): {t_peak:.2f} days")
        
        # Maximum concentration C_max in g/m^3
        # Term 1: Mass loading and dispersion
        term1 = MTotal / math.sqrt(4 * math.pi * D * t_peak)
        # Term 2: Advection-dispersion exponent
        exp_numerator = -((x - v * t_peak)**2)
        exp_denominator = 4 * D * t_peak
        term2 = math.exp(exp_numerator / exp_denominator)
        # Term 3: Decay
        term3 = math.exp(-k * t_peak)
        # Full equation for C_max in g/m^3
        C_max_g_m3 = term1 * term2 * term3 * F
        
        # Convert to ug/L (1 g/m^3 = 1000 ug/L)
        C_max_ug_L = C_max_g_m3 * 1000
        print(f"  - Peak Concentration (C_max): {C_max_ug_L:.3e} ug/L")
        
        # Hazard Quotient HQ
        hq = C_max_ug_L / props['ec50_ug_L']
        print(f"  - Freshwater Algal EC50: {props['ec50_ug_L']} ug/L")
        print(f"  - Hazard Quotient (HQ): {hq:.3e}")
        
        results.append({'name': name, 'C_max': C_max_ug_L, 'HQ': hq})

    # --- 3. Process and Print Final Results ---
    # Find the chemical with the highest concentration
    highest_conc_chem = max(results, key=lambda r: r['C_max'])
    
    print("\n" + "="*50)
    print("                      FINAL RESULTS")
    print("="*50)

    print(f"\nThe highest concentration reaching the spring is from {highest_conc_chem['name']}.")
    print(f"Highest Individual Concentration: {highest_conc_chem['C_max']:.4e} ug/L")

    # Calculate and discuss mixture effect using Hazard Index (HI)
    hazard_index = sum(res['HQ'] for res in results)
    
    print("\n--- Mixture Effect Analysis (assuming additivity) ---")
    print(f"The Hazard Index (HI) is calculated as the sum of Hazard Quotients (HQ = C_max / EC50).")
    for res in results:
        print(f"  HQ for {res['name']}: {res['C_max']:.2e} / {chemicals_data[res['name']]['ec50_ug_L']} = {res['HQ']:.2e}")
    
    print(f"\nTotal Hazard Index (HI) = {results[0]['HQ']:.2e} (Atrazine) + {results[1]['HQ']:.2e} (PFOS) + {results[2]['HQ']:.2e} (Endosulfan)")
    print(f"Total Hazard Index (HI) = {hazard_index:.4e}")

    if hazard_index < 1:
        print("\nEffect of the Mixture: ADDITIVE (by assumption), but NEGLIGIBLE.")
        print("The Hazard Index is significantly less than 1. This indicates that even when assuming the toxic effects are additive,")
        print("the combined concentration of the chemicals reaching the spring is far below a level expected to cause harm to the algal community.")
    else:
        print("\nEffect of the Mixture: Additive and Potentially Significant.")
        print("The Hazard Index is greater than or equal to 1, indicating a potential risk to the algal community from the additive effect of the mixture.")


# Run the main function
solve_transport_problem()
<<<1.8353e-05>>>