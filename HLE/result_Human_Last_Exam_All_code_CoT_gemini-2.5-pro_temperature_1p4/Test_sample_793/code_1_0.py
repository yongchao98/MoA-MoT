import math

def solve_environmental_transport():
    # --- 1. Initial Parameters ---
    # Scenario parameters
    V_water_L = 1500000
    rho_coal_kg_m3 = 1346
    f_oc = 0.50
    theta = 0.30
    fracture_length_m = 80
    fracture_height_m = 10
    distance_x_m = 100
    advection_v_md = 0.5
    dispersivity_alpha_m = 0.5

    # Chemical properties dictionary
    chemicals = {
        "Atrazine": {
            "prod_percent": 0.01, "conc_in_prod_ug_L": 40, "log_koc": 2.20,
            "half_life_d": 90, "ec50_ug_L": 100
        },
        "PFOS": {
            "prod_percent": 0.001, "conc_in_prod_ug_L": 300, "log_koc": 3.65,
            "half_life_d": 14965, "ec50_ug_L": 480
        },
        "Endosulfan": {
            "prod_percent": 0.005, "conc_in_prod_ug_L": 20, "log_koc": 4.3,
            "half_life_d": 60, "ec50_ug_L": 560
        },
    }

    # --- 2. Calculate Shared Parameters ---
    rho_bulk_kg_m3 = rho_coal_kg_m3 * (1 - theta)
    dispersion_D_m2d = dispersivity_alpha_m * advection_v_md
    fracture_area_m2 = fracture_length_m * fracture_height_m

    print("--- Shared Aquifer and Transport Parameters ---")
    print(f"Coal Bulk Density (ρ_b): {rho_bulk_kg_m3:.2f} kg/m³")
    print(f"Dispersion Coefficient (D): {dispersion_D_m2d} m²/d")
    print(f"Fracture Face Area (A): {fracture_area_m2} m²")
    print(f"Porosity (θ): {theta}")
    print("-" * 50 + "\n")

    results = {}
    hazard_index = 0

    # --- 3. Process Each Chemical ---
    for name, props in chemicals.items():
        print(f"--- Calculations for {name} ---")

        # Mass calculation
        m_total_ug = V_water_L * props["prod_percent"] * props["conc_in_prod_ug_L"]
        print(f"Total Mass (M_total): {m_total_ug:,.0f} µg")

        # Partitioning and Retardation
        koc_L_kg = 10**props["log_koc"]
        kd_L_kg = koc_L_kg * f_oc
        kd_m3_kg = kd_L_kg / 1000  # Convert L/kg to m³/kg
        R = 1 + (rho_bulk_kg_m3 / theta) * kd_m3_kg
        print(f"Partition Coefficient (Kd): {kd_L_kg:.2f} L/kg")
        print(f"Retardation Factor (R): {R:.2f}")

        # Decay
        k_per_d = math.log(2) / props["half_life_d"]
        print(f"Decay Rate (k): {k_per_d:.6f} /day")

        # Peak Time (t_peak) calculation by solving quadratic equation At^2 + Bt + C = 0
        a = advection_v_md**2 + 4 * dispersion_D_m2d * k_per_d
        b = 2 * dispersion_D_m2d
        c = -distance_x_m**2
        t_peak_d = (-b + math.sqrt(b**2 - 4 * a * c)) / (2 * a)
        print(f"Peak Travel Time (t_peak): {t_peak_d:.2f} days")

        # Peak Concentration (C_max) calculation
        # Formula: C(x,t) = [M_total / (A * θ * sqrt(4πDt))] * exp(-(x-vt)²/(4Dt)) * exp(-kt) * (1/R)
        
        # Components of the equation
        sqrt_4_pi_Dt = math.sqrt(4 * math.pi * dispersion_D_m2d * t_peak_d)
        exp_advection_term_exponent = -((distance_x_m - advection_v_md * t_peak_d)**2) / (4 * dispersion_D_m2d * t_peak_d)
        exp_advection_term = math.exp(exp_advection_term_exponent)
        exp_decay_term = math.exp(-k_per_d * t_peak_d)
        inv_R_term = 1 / R

        # Final concentration in µg/m³
        c_max_ug_m3 = (m_total_ug / (fracture_area_m2 * theta * sqrt_4_pi_Dt)) * exp_advection_term * exp_decay_term * inv_R_term
        # Convert to µg/L
        c_max_ug_L = c_max_ug_m3 / 1000

        print(f"Peak Concentration (C_max): {c_max_ug_L:.3e} µg/L")

        # Hazard Quotient
        hq = c_max_ug_L / props["ec50_ug_L"]
        hazard_index += hq
        print(f"Hazard Quotient (HQ = C_max/EC50): {hq:.3e}")
        
        results[name] = {'c_max_ug_L': c_max_ug_L, 't_peak_d': t_peak_d, 'k_per_d': k_per_d, 'R': R, 'm_total_ug': m_total_ug}
        print("-" * 50 + "\n")

    # --- 4. Identify Highest Concentration ---
    highest_conc_chem = max(results, key=lambda k: results[k]['c_max_ug_L'])
    highest_conc_val = results[highest_conc_chem]['c_max_ug_L']

    print("--- Final Results ---")
    print(f"Highest concentration at the spring is from: {highest_conc_chem}")
    print(f"Highest concentration value: {highest_conc_val:.3e} µg/L")
    
    # --- 5. Assess Mixture Effect ---
    print(f"\nMixture Effect Assessment:")
    print("The effect of the mixture is assessed by summing the Hazard Quotients (HQ) of each chemical, which assumes an additive effect.")
    print(f"Total Hazard Index (HI = Σ HQ): {hazard_index:.3e}")
    if hazard_index < 1:
        print("Since the Hazard Index is much less than 1, the additive effect is considered negligible.")
    else:
        print("Since the Hazard Index is 1 or greater, the additive effect may be of concern.")

    # --- 6. Print Final Equation with Numbers ---
    print("\n--- Final Equation for Highest Concentration Chemical ({}) ---".format(highest_conc_chem))
    chem_res = results[highest_conc_chem]
    print("C_max = [M_total / (A * θ * sqrt(4*π*D*t))] * exp(-(x-v*t)²/(4*D*t)) * exp(-k*t) * (1/R)")
    print("C_max = [{m_total:,.0f} µg / ({A} m² * {theta} * sqrt(4*π*{D} m²/d * {t:.2f} d))] * exp(-({x} m - {v} m/d * {t:.2f} d)² / (4*{D} m²/d * {t:.2f} d)) * exp(-{k:.6f} /d * {t:.2f} d) * (1/{R:.2f})".format(
        m_total=chem_res['m_total_ug'], A=fracture_area_m2, theta=theta, D=dispersion_D_m2d, t=chem_res['t_peak_d'],
        x=distance_x_m, v=advection_v_md, k=chem_res['k_per_d'], R=chem_res['R']
    ))
    
    final_answer = (f"The highest concentration reaching the spring is {highest_conc_val:.3e} µg/L from {highest_conc_chem}. "
                    f"The effect of the mixture is additive. The calculated Hazard Index is {hazard_index:.3e}, "
                    "which is far below the threshold of 1, indicating a negligible risk to the algal community.")

    print(f"\n<<<{final_answer}>>>")

# Run the simulation
solve_environmental_transport()
