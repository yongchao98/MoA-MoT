import numpy as np

def solve_contamination_scenario():
    """
    Calculates contaminant transport and risk from a hydraulic fracturing operation.
    """

    # --- 1. Define Inputs ---

    # Environmental Parameters
    V_injected_total_L = 1_500_000  # L
    coal_density_kg_m3 = 1346  # kg/m^3
    fracture_organic_carbon = 0.50  # 50%
    fracture_porosity = 0.30  # 30%
    fracture_L = 80  # m
    fracture_H = 10  # m
    fracture_W = 0.005 # m
    distance_to_spring_m = 100  # m
    advection_velocity_m_d = 0.5  # m/d
    dispersivity_m = 0.5 # m

    # Chemical properties in a list of dictionaries
    chemicals = [
        {
            'name': 'Atrazine',
            'product_percent': 0.01, # 1%
            'conc_in_prod_ug_L': 40,
            'log_koc': 2.20,
            'half_life_d': 90,
            'ec50_ug_L': 100
        },
        {
            'name': 'PFOS',
            'product_percent': 0.001, # 0.1%
            'conc_in_prod_ug_L': 300,
            'log_koc': 3.65,
            'half_life_d': 14965,
            'ec50_ug_L': 480
        },
        {
            'name': 'Endosulfan',
            'product_percent': 0.005, # 0.5%
            'conc_in_prod_ug_L': 20,
            'log_koc': 4.30,
            'half_life_d': 60,
            'ec50_ug_L': 560
        }
    ]

    # --- 2. Calculate Intermediate and Model Parameters ---
    
    # Convert density for consistent units (kg/m^3 -> kg/L)
    bulk_density_kg_L = coal_density_kg_m3 / 1000

    # Longitudinal Dispersion Coefficient (D_L)
    D_L = dispersivity_m * advection_velocity_m_d  # m^2/d

    # Area of the source plane (fracture face)
    A_source = fracture_L * fracture_H # m^2

    # A list to store detailed results for each chemical
    results = []

    for chem in chemicals:
        # Mass of contaminant in micrograms (ug)
        volume_product_L = V_injected_total_L * chem['product_percent']
        mass_ug = volume_product_L * chem['conc_in_prod_ug_L']
        chem['mass_ug'] = mass_ug
        
        # Partitioning coefficients
        koc = 10**chem['log_koc']  # L/kg
        kd = koc * fracture_organic_carbon  # L/kg
        
        # Retardation Factor (R)
        R = 1 + (bulk_density_kg_L / fracture_porosity) * kd
        chem['R'] = R

        # First-order decay rate (k)
        k = np.log(2) / chem['half_life_d']  # 1/day
        chem['k'] = k
        
        results.append(chem)

    # --- 3. Concentration Calculation Functions and Peak Finding ---
    
    def calculate_concentration(t, chem_params):
        """Calculates concentration C(x,t) based on the chosen model."""
        if t <= 0:
            return 0
        
        M_A = chem_params['mass_ug'] / A_source # Mass per area (ug/m^2)
        n = fracture_porosity
        x = distance_to_spring_m
        v = advection_velocity_m_d
        D = D_L
        k = chem_params['k']
        R = chem_params['R']

        # Numerically stable calculation
        sqrt_term = np.sqrt(4 * np.pi * D * t)
        
        # Handle potential division by zero for very small t
        if sqrt_term == 0:
            return 0

        exp1_arg = -((x - v * t)**2) / (4 * D * t)
        exp2_arg = -k * t

        # Avoid overflow for very large negative exponents
        if exp1_arg < -700 or exp2_arg < -700:
             return 0

        conc_ug_m3 = (M_A / n) / sqrt_term * np.exp(exp1_arg) * np.exp(exp2_arg) * (1 / R)
        
        # Convert from ug/m^3 to ug/L
        conc_ug_L = conc_ug_m3 / 1000
        return conc_ug_L

    # Find the peak concentration for each chemical
    for chem in results:
        # To find the time of the peak, solve the quadratic equation:
        # t^2*(v^2 + 4Dk) + t*(2D) - x^2 = 0
        a = advection_velocity_m_d**2 + 4 * D_L * chem['k']
        b = 2 * D_L
        c = -distance_to_spring_m**2
        
        # Using the quadratic formula for t > 0
        t_peak = (-b + np.sqrt(b**2 - 4 * a * c)) / (2 * a)
        
        chem['t_peak_d'] = t_peak
        chem['c_max_ug_L'] = calculate_concentration(t_peak, chem)
        
    # --- 4. Identify Highest Concentration and Time of Event ---
    
    primary_contaminant = max(results, key=lambda x: x['c_max_ug_L'])
    highest_conc = primary_contaminant['c_max_ug_L']
    time_of_event = primary_contaminant['t_peak_d']

    print(f"The highest concentration of an individual chemical reaching the spring is {highest_conc:.6f} μg/L.")
    print(f"This peak concentration is from {primary_contaminant['name']} and occurs after {time_of_event:.1f} days.")
    print("-" * 30)

    # --- 5. Assess Mixture Effect at Time of Event ---

    print(f"Assessing the mixture effect at t = {time_of_event:.1f} days:")
    
    hazard_index = 0
    concentrations_at_event_time = {}

    for chem in results:
        conc_at_event = calculate_concentration(time_of_event, chem)
        concentrations_at_event_time[chem['name']] = conc_at_event
        hazard_index += conc_at_event / chem['ec50_ug_L']

    print("Concentrations of each chemical at this time:")
    for name, conc in concentrations_at_event_time.items():
        print(f"- {name}: {conc:.6f} μg/L")
    
    print("-" * 30)
    print("The effect of the mixture is assessed using the Hazard Index (HI).")
    print("HI is the sum of each chemical's concentration divided by its EC50 value.")
    
    print("\nHazard Index (HI) = (C_Atrazine / EC50_Atrazine) + (C_PFOS / EC50_PFOS) + (C_Endosulfan / EC50_Endosulfan)")
    
    c_a = concentrations_at_event_time['Atrazine']
    c_p = concentrations_at_event_time['PFOS']
    c_e = concentrations_at_event_time['Endosulfan']
    
    ec50_a = next(c for c in results if c['name'] == 'Atrazine')['ec50_ug_L']
    ec50_p = next(c for c in results if c['name'] == 'PFOS')['ec50_ug_L']
    ec50_e = next(c for c in results if c['name'] == 'Endosulfan')['ec50_ug_L']
    
    print(f"HI = ({c_a:.6f} / {ec50_a}) + ({c_p:.6f} / {ec50_p}) + ({c_e:.6f} / {ec50_e})")
    print(f"HI = {c_a/ec50_a:.6e} + {c_p/ec50_p:.6e} + {c_e/ec50_e:.6e}")
    print(f"HI = {hazard_index:.6f}")

    if hazard_index >= 1:
        effect = "additive and likely to pose a risk to the algal community."
    else:
        effect = "additive, but unlikely to pose a significant risk to the algal community at this concentration."

    print(f"\nBased on a Hazard Index of {hazard_index:.6f}, the combined effect of the mixture is considered {effect}")


solve_contamination_scenario()