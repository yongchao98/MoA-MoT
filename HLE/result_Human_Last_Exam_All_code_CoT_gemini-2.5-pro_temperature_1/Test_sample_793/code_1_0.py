import numpy as np

def solve_contaminant_transport():
    """
    Calculates the peak contaminant concentration at a spring and assesses mixture effects
    using the Advection-Dispersion Equation.
    """
    # --- 1. Define constants from the problem statement ---
    # Scenario constants
    V_water_total_L = 1_500_000
    rho_coal_particle_kg_m3 = 1346
    f_oc = 0.50  # Organic carbon fraction
    theta = 0.30  # Porosity / water content
    x_dist_m = 100
    v_gw_m_d = 0.5  # Advection velocity
    dispersivity_m = 0.5
    
    # The flow is 1D away from the fracture face. The cross-sectional area (A)
    # is the area of this face.
    fracture_area_A_m2 = 80 * 10

    # Chemical-specific data
    chemicals = {
        'Atrazine': {
            'prod_frac': 0.01, 'conc_in_prod_ug_L': 40, 'log_koc': 2.20,
            'half_life_d': 90, 'ec50_ug_L': 100
        },
        'PFOS': {
            'prod_frac': 0.001, 'conc_in_prod_ug_L': 300, 'log_koc': 3.65,
            'half_life_d': 14965, 'ec50_ug_L': 480
        },
        'Endosulfan': {
            'prod_frac': 0.005, 'conc_in_prod_ug_L': 20, 'log_koc': 4.30,
            'half_life_d': 60, 'ec50_ug_L': 560
        }
    }

    # --- 2. Calculate derived parameters ---
    # Bulk density (kg/L)
    rho_b_kg_m3 = rho_coal_particle_kg_m3 * (1 - theta)
    rho_b_kg_L = rho_b_kg_m3 / 1000
    # Dispersion coefficient (m^2/d)
    D_disp_m2_d = dispersivity_m * v_gw_m_d

    # Calculate parameters for each chemical
    for name, params in chemicals.items():
        # Total mass injected (ug)
        params['MTotal_ug'] = (V_water_total_L * params['prod_frac']) * params['conc_in_prod_ug_L']
        # Koc (L/kg) and Kd (L/kg)
        params['Koc_L_kg'] = 10**params['log_koc']
        params['Kd_L_kg'] = params['Koc_L_kg'] * f_oc
        # Decay rate k (1/d)
        params['k_per_d'] = np.log(2) / params['half_life_d']
        # Retardation factor R (dimensionless)
        params['R'] = 1 + (rho_b_kg_L * params['Kd_L_kg']) / theta

    # --- 3. Define the Advection-Dispersion Equation (ADE) ---
    def advection_dispersion_1d(t, M, R, D, v, k, A, n, x):
        if t <= 0: return 0
        
        # Account for retardation on velocity and dispersion
        v_retarded = v / R
        D_retarded = D / R
        
        # Avoid division by zero at t=0
        if D_retarded * t == 0: return 0.0

        # Pre-exponential factor
        term1_denom = A * n * np.sqrt(4 * np.pi * D_retarded * t)
        if term1_denom == 0: return 0.0
        term1 = M / term1_denom
        
        # Exponential term for advection and dispersion
        exp_term_numerator = -(x - v_retarded * t)**2
        exp_term_denominator = 4 * D_retarded * t
        exp_term = np.exp(exp_term_numerator / exp_term_denominator)
        
        # Exponential term for decay
        decay_term = np.exp(-k * t)
        
        # Final concentration in ug/m^3
        conc_ug_m3 = term1 * exp_term * decay_term
        return conc_ug_m3

    # --- 4. Find the peak concentration for each chemical by simulating over time ---
    # Time range needs to be large enough to capture highly retarded peaks
    time_points_d = np.logspace(3, 7, 2000) # Simulate from 1,000 to 10,000,000 days
    
    peak_results = {}
    highest_overall_conc_ug_L = -1
    chem_with_highest_conc = None
    time_of_highest_conc_d = -1

    for name, params in chemicals.items():
        # Calculate concentration over time
        concentrations_ug_m3 = [advection_dispersion_1d(
            t=t, M=params['MTotal_ug'], R=params['R'], D=D_disp_m2_d, v=v_gw_m_d,
            k=params['k_per_d'], A=fracture_area_A_m2, n=theta, x=x_dist_m
        ) for t in time_points_d]
        
        # Find the peak concentration and time for this chemical
        peak_conc_ug_m3 = np.max(concentrations_ug_m3)
        peak_conc_ug_L = peak_conc_ug_m3 / 1000 # Convert m^3 to L
        
        # Check if this is the highest peak found so far
        if peak_conc_ug_L > highest_overall_conc_ug_L:
            highest_overall_conc_ug_L = peak_conc_ug_L
            chem_with_highest_conc = name
            time_of_highest_conc_d = time_points_d[np.argmax(concentrations_ug_m3)]

    print(f"The highest concentration of an individual chemical reaching the spring is from '{chem_with_highest_conc}'.")
    print(f"Highest Concentration: {highest_overall_conc_ug_L:.6f} μg/L")
    print(f"This peak is predicted to occur at approximately {time_of_highest_conc_d:.0f} days.")
    print("-" * 50)
    
    # To satisfy the prompt's requirement, here is the final equation with numbers for the highest concentration found
    p = chemicals[chem_with_highest_conc]
    v_r = v_gw_m_d / p['R']
    D_r = D_disp_m2_d / p['R']
    print("Calculation for the highest concentration:")
    print(f"C(x={x_dist_m}, t={time_of_highest_conc_d:.0f}) = [M={p['MTotal_ug']:.2e} / (A={fracture_area_A_m2} * n={theta} * R={p['R']:.2f})] * "
          f"[1 / sqrt(4 * pi * D_r={D_r:.4f} * t={time_of_highest_conc_d:.0f})] * "
          f"exp[-((x={x_dist_m} - v_r={v_r:.4f} * t={time_of_highest_conc_d:.0f})^2 / "
          f"(4 * D_r={D_r:.4f} * t={time_of_highest_conc_d:.0f}))] * exp[-(k={p['k_per_d']:.4f} * t={time_of_highest_conc_d:.0f})]")
    print("-" * 50)

    # --- 5. Analyze mixture effect at the time of the highest peak ---
    print(f"Mixture effect analysis at t = {time_of_highest_conc_d:.0f} days:")
    sum_of_rq = 0
    for name, params in chemicals.items():
        # Calculate concentration of each chemical at that specific time
        conc_at_peak_time_ug_m3 = advection_dispersion_1d(
            t=time_of_highest_conc_d, M=params['MTotal_ug'], R=params['R'], D=D_disp_m2_d,
            v=v_gw_m_d, k=params['k_per_d'], A=fracture_area_A_m2, n=theta, x=x_dist_m
        )
        conc_at_peak_time_ug_L = conc_at_peak_time_ug_m3 / 1000
        
        # Calculate Risk Quotient (RQ)
        rq = conc_at_peak_time_ug_L / params['ec50_ug_L']
        sum_of_rq += rq
        print(f"  {name}: {conc_at_peak_time_ug_L:.6f} μg/L (RQ = {rq:.4f})")

    print(f"\nSum of Risk Quotients (RQ_mixture) = {sum_of_rq:.4f}")
    
    # The default assumption for mixture toxicity is additivity.
    effect = "additive"
    if sum_of_rq > 1:
        print(f"Conclusion: The mixture poses a potential risk based on an {effect} effect.")
    else:
        print(f"Conclusion: The mixture is unlikely to pose a risk based on an {effect} effect.")
        
    # Format the final answer as requested
    final_answer = f"{highest_overall_conc_ug_L:.6f}, {effect}"
    print(f"\n<<<{final_answer}>>>")

solve_contaminant_transport()