import numpy as np

def solve_transport_problem():
    """
    Calculates the highest chemical concentration reaching a spring and assesses
    the mixture toxicity based on the provided environmental parameters.
    """
    # --- 1. Define Constants and Inputs ---

    # General parameters
    V_total_water_L = 1_500_000  # L
    rho_b_kg_m3 = 1346      # Bulk density of coal seam (kg/m^3)
    f_oc = 0.50             # Fraction of organic carbon in coal seam
    theta = 0.30            # Porosity of coal seam
    x_m = 100               # Distance to spring (m)
    v_gw_m_d = 0.5          # Groundwater advection velocity (m/d)
    alpha_m = 0.5           # Dispersivity (m)
    fracture_h_m = 10       # Fracture height (m)
    fracture_l_m = 80       # Fracture length (m)
    
    # Cross-sectional area of flow, assumed to be the fracture face
    A_m2 = fracture_h_m * fracture_l_m 

    # Convert bulk density from kg/m^3 to kg/L for retardation calculation
    rho_b_kg_L = rho_b_kg_m3 / 1000.0

    # Hydrodynamic dispersion coefficient (D = alpha * v)
    D_m2_d = alpha_m * v_gw_m_d

    # Chemical-specific data
    chemicals = {
        "Atrazine": {
            "product_fraction": 0.01, "conc_ug_L": 40, "log_koc": 2.20,
            "half_life_d": 90, "ec50_ug_L": 100,
        },
        "PFOS": {
            "product_fraction": 0.001, "conc_ug_L": 300, "log_koc": 3.65,
            "half_life_d": 14965, "ec50_ug_L": 480,
        },
        "Endosulfan": {
            "product_fraction": 0.005, "conc_ug_L": 20, "log_koc": 4.3,
            "half_life_d": 60, "ec50_ug_L": 560,
        }
    }

    # --- 2. Calculate Initial Mass and Transport Parameters ---
    for name, props in chemicals.items():
        # Initial Mass (M) in micrograms
        V_product_L = V_total_water_L * props["product_fraction"]
        props["M_ug"] = V_product_L * props["conc_ug_L"]
        
        # Koc and Kd
        props["koc_L_kg"] = 10**props["log_koc"]
        props["kd_L_kg"] = props["koc_L_kg"] * f_oc
        
        # Retardation Factor (R)
        props["R"] = 1 + (rho_b_kg_L / theta) * props["kd_L_kg"]
        
        # Retarded velocity (v_R) and dispersion (D_R)
        props["v_R_m_d"] = v_gw_m_d / props["R"]
        props["D_R_m2_d"] = D_m2_d / props["R"]
        
        # First-order decay rate (k)
        props["k_per_d"] = np.log(2) / props["half_life_d"]

    # --- 3. Define Concentration Function and Find Peak ---

    def get_concentration_at_time(t_d, chem_props):
        """Calculates concentration C(x,t) in ug/L for a given chemical."""
        if t_d <= 1e-9: return 0.0 # Avoid division by zero at t=0

        M = chem_props["M_ug"]
        v_R = chem_props["v_R_m_d"]
        D_R = chem_props["D_R_m2_d"]
        k = chem_props["k_per_d"]
        
        # Standard 1D Advection-Dispersion-Reaction Equation
        term1 = M / (A_m2 * theta)
        term2 = 1.0 / np.sqrt(4 * np.pi * D_R * t_d)
        exponent_val = -((x_m - v_R * t_d)**2) / (4 * D_R * t_d)
        term3 = np.exp(exponent_val)
        term4 = np.exp(-k * t_d)
        
        conc_ug_m3 = term1 * term2 * term3 * term4
        return conc_ug_m3 / 1000.0 # Convert from ug/m^3 to ug/L

    peak_concentrations = {}
    peak_times = {}

    for name, props in chemicals.items():
        # Estimate peak time to create a search window
        t_peak_approx = x_m / props["v_R_m_d"]
        
        # Create a time array to search for the peak concentration
        time_array = np.linspace(t_peak_approx * 0.2, t_peak_approx * 1.8, 5000)
        
        concentrations = np.array([get_concentration_at_time(t, props) for t in time_array])
        
        max_conc = np.max(concentrations)
        max_time = time_array[np.argmax(concentrations)]
        
        peak_concentrations[name] = max_conc
        peak_times[name] = max_time

    # --- 4. Determine Highest Concentration and Mixture Effect ---

    highest_conc_chem = max(peak_concentrations, key=peak_concentrations.get)
    highest_conc_value = peak_concentrations[highest_conc_chem]
    
    # The standard assumption for mixture toxicity without further data is additivity.
    mixture_effect = "additive"

    # --- 5. Print Final Results ---
    
    print("="*50)
    print("Highest Concentration Analysis")
    print("="*50)
    print(f"The chemical predicted to have the highest concentration at the spring is: {highest_conc_chem}")
    
    # Output the numbers used in the final equation for this chemical
    props = chemicals[highest_conc_chem]
    t_peak = peak_times[highest_conc_chem]
    
    print("\nThe calculation for the peak concentration C(x,t) uses the following parameters:")
    print(f"C(t) = (M / (A * n)) * (1/sqrt(4*pi*D_R*t)) * exp(-(x-v_R*t)^2/(4*D_R*t)) * exp(-k*t) / 1000")
    print("\n--- Equation Parameters ---")
    print(f"M (mass)              = {props['M_ug']:.3e} ug")
    print(f"A (area)                = {A_m2:.1f} m^2")
    print(f"n (porosity)            = {theta:.2f}")
    print(f"D_R (retarded dispersion) = {props['D_R_m2_d']:.3e} m^2/d")
    print(f"v_R (retarded velocity)   = {props['v_R_m_d']:.3e} m/d")
    print(f"k (decay rate)          = {props['k_per_d']:.3e} 1/d")
    print(f"x (distance)            = {x_m:.1f} m")
    print(f"t (time of peak)        = {t_peak:.2e} days")
    print("---------------------------\n")

    print("Plugging the values into the equation at the time of the peak:")
    print(f"C_peak = ({props['M_ug']:.3e} / ({A_m2:.1f} * {theta:.2f})) * "
          f"(1/sqrt(4*pi*{props['D_R_m2_d']:.3e}*{t_peak:.2e})) * "
          f"exp(-({x_m:.1f}-{props['v_R_m_d']:.3e}*{t_peak:.2e})^2/(4*{props['D_R_m2_d']:.3e}*{t_peak:.2e})) * "
          f"exp(-{props['k_per_d']:.3e}*{t_peak:.2e}) / 1000")

    print("\n--- Final Answer ---")
    print(f"Highest concentration of an individual chemical: {highest_conc_value:.4e} ug/L (from {highest_conc_chem})")
    print(f"Effect of the mixture at this concentration: {mixture_effect}")
    print("="*50)

# Run the analysis
solve_transport_problem()