import math

def solve_contaminant_transport():
    # Step 1: Define Inputs and Constants
    # Using consistent units: meters, days, micrograms (µg)

    # General Parameters
    v_total_water_litres = 1500000
    rho_coal_particle = 1346  # kg/m³
    f_oc = 0.50  # fraction of organic carbon
    porosity_n = 0.30  # fraction of water
    distance_x = 100  # m
    advection_velocity_v = 0.5  # m/d
    dispersivity_alpha_L = 0.5  # m
    fracture_length = 80 # m
    fracture_height = 10 # m

    # Derived constants
    cross_sectional_area_A = fracture_length * fracture_height # m²
    dispersion_coeff_D = dispersivity_alpha_L * advection_velocity_v # m²/d
    bulk_density_rhob = (1 - porosity_n) * rho_coal_particle  # kg/m³
    
    # List of dictionaries to hold chemical-specific data
    chemicals = [
        {
            'name': 'Atrazine',
            'product_fraction': 0.01,
            'concentration_ug_L': 40,
            'log_koc': 2.20,
            'half_life_days': 90,
            'ec50_ug_L': 100
        },
        {
            'name': 'PFOS',
            'product_fraction': 0.001,
            'concentration_ug_L': 300,
            'log_koc': 3.65,
            'half_life_days': 14965,
            'ec50_ug_L': 480
        },
        {
            'name': 'Endosulfan',
            'product_fraction': 0.005,
            'concentration_ug_L': 20,
            'log_koc': 4.3,
            'half_life_days': 60,
            'ec50_ug_L': 560
        }
    ]

    results = []
    hazard_quotients = []

    # Step 2 & 3: Calculate parameters and transport for each chemical
    for chem in chemicals:
        # Calculate Total Mass (M) in µg
        # M = total water (L) * product fraction * concentration (µg/L)
        mass_M = v_total_water_litres * chem['product_fraction'] * chem['concentration_ug_L']
        chem['mass_M'] = mass_M
        
        # Mass loading per unit area (M/A) in µg/m²
        mass_loading_M_A = mass_M / cross_sectional_area_A
        chem['mass_loading_M_A'] = mass_loading_M_A

        # Calculate first-order decay rate constant (k) in 1/d
        k = math.log(2) / chem['half_life_days']
        chem['k'] = k
        
        # Calculate Koc and Kd
        koc = 10**chem['log_koc']  # L/kg
        # Convert Kd to m³/kg for consistency with rho_b (kg/m³)
        kd_m3_kg = koc * f_oc * 0.001 # m³/kg
        chem['kd_m3_kg'] = kd_m3_kg

        # Calculate Retardation Factor (R)
        R = 1 + (bulk_density_rhob * kd_m3_kg) / porosity_n
        chem['R'] = R

        # Step 4: Solve for Peak Concentration using the Simplified Model
        # Find time to peak for unretarded case (R=1) by solving quadratic equation
        # for t: t^2*(v^2 + 4Dk) + t*(2D) - x^2 = 0
        a_quad = advection_velocity_v**2 + 4 * dispersion_coeff_D * k
        b_quad = 2 * dispersion_coeff_D
        c_quad = -distance_x**2
        t_peak_unretarded = (-b_quad + math.sqrt(b_quad**2 - 4 * a_quad * c_quad)) / (2 * a_quad)
        chem['t_peak_unretarded'] = t_peak_unretarded
        
        # Calculate peak concentration for unretarded case C_peak(R=1)
        # C(x,t) = (M/A) / (n * sqrt(4*pi*D*t)) * exp(-(x-vt)^2/(4Dt)) * exp(-kt)
        term1_denom = porosity_n * math.sqrt(4 * math.pi * dispersion_coeff_D * t_peak_unretarded)
        term1 = mass_loading_M_A / term1_denom

        exp_arg1_num = -(distance_x - advection_velocity_v * t_peak_unretarded)**2
        exp_arg1_den = 4 * dispersion_coeff_D * t_peak_unretarded
        exp_arg2 = -k * t_peak_unretarded

        C_peak_unretarded_ug_m3 = term1 * math.exp(exp_arg1_num / exp_arg1_den) * math.exp(exp_arg2)
        chem['C_peak_unretarded_ug_m3'] = C_peak_unretarded_ug_m3

        # Apply retardation factor and convert units from µg/m³ to µg/L
        C_peak_final_ug_m3 = C_peak_unretarded_ug_m3 / R
        C_peak_final_ug_L = C_peak_final_ug_m3 / 1000
        chem['C_peak_final_ug_L'] = C_peak_final_ug_L
        results.append(chem)
        
        # Step 5: Assess Mixture Effect
        hq = C_peak_final_ug_L / chem['ec50_ug_L']
        chem['hq'] = hq
        hazard_quotients.append(hq)

    # Find the chemical with the highest concentration
    highest_conc_chem = max(results, key=lambda x: x['C_peak_final_ug_L'])
    
    # Calculate the Hazard Index (HI)
    hazard_index = sum(hazard_quotients)
    
    # Step 6: Final Output
    print("--- Analysis of Contaminant Transport ---")
    print(f"The highest concentration reaching the spring belongs to {highest_conc_chem['name']}.")
    print("\nDetailed calculation for this chemical:")
    
    # Print the equation breakdown
    print(f"\n1. Unretarded Peak Concentration C_peak(R=1) in µg/m³:")
    print("   C_peak(R=1) = (M/A) / (n * sqrt(4*pi*D*t)) * exp(-(x-vt)²/(4Dt)) * exp(-kt)")
    M_A_val = highest_conc_chem['mass_loading_M_A']
    n_val = porosity_n
    D_val = dispersion_coeff_D
    t_val = highest_conc_chem['t_peak_unretarded']
    x_val = distance_x
    v_val = advection_velocity_v
    k_val = highest_conc_chem['k']
    print(f"   C_peak(R=1) = ({M_A_val:.2f} / ({n_val} * sqrt(4 * pi * {D_val:.2f} * {t_val:.2f}))) * exp(-({x_val} - {v_val}*{t_val:.2f})² / (4*{D_val:.2f}*{t_val:.2f})) * exp(-{k_val:.4f}*{t_val:.2f})")
    print(f"   C_peak(R=1) = {highest_conc_chem['C_peak_unretarded_ug_m3']:.4f} µg/m³")

    print(f"\n2. Final Peak Concentration (applying retardation R = {highest_conc_chem['R']:.2f}):")
    print("   C_final = C_peak(R=1) / R")
    C_unretarded_val = highest_conc_chem['C_peak_unretarded_ug_m3']
    R_val = highest_conc_chem['R']
    print(f"   C_final = {C_unretarded_val:.4f} / {R_val:.2f} = {C_unretarded_val / R_val:.4f} µg/m³")
    
    final_conc_ug_L = highest_conc_chem['C_peak_final_ug_L']
    print(f"   C_final = {final_conc_ug_L * 1000:.4f} µg/m³ * (1 m³ / 1000 L) = {final_conc_ug_L:.4e} µg/L")
    
    print("\n--- Final Answer ---")
    print(f"\nHighest Concentration: The highest concentration of an individual chemical reaching the spring is {final_conc_ug_L:.4e} µg/L, which is for {highest_conc_chem['name']}.")
    
    # Assess mixture effect based on Hazard Index
    effect = "additive"
    print(f"\nMixture Effect: Based on a Hazard Index (HI) of {hazard_index:.4e}, the potential for harm from the chemical mixture is exceedingly low.")
    print(f"In the absence of specific interaction data, the mixture effect is assumed to be {effect}.")


solve_contaminant_transport()
<<<8.7238e-05>>>