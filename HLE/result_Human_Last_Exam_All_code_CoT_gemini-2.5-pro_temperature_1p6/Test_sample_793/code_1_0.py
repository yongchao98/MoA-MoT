import math
import scipy.optimize as optimize

def solve_transport_problem():
    """
    Calculates the peak concentration of contaminants reaching a spring
    and assesses the combined toxic effect.
    """
    # 1. Define constants and input parameters
    # Environmental parameters
    total_water_L = 1500000
    coal_density_kg_m3 = 1346  # Solid particle density
    coal_water_fraction_mass = 0.3
    coal_foc = 0.5 # fraction of organic carbon
    distance_to_spring_m = 100
    fracture_length_m = 80
    fracture_height_m = 10
    advection_velocity_m_d = 0.5
    dispersivity_m = 0.5

    # Chemical properties dictionary
    # M_product_%: Percentage of total injected water volume
    # M_conc_ug_L: Concentration in the product
    # log_koc: Log of organic carbon-water partitioning coefficient
    # t_half_life_d: Half-life in days
    # ec50_ug_L: Freshwater algal EC50
    chemicals_data = {
        'Atrazine': {
            'M_product_%': 0.01,
            'M_conc_ug_L': 40,
            'log_koc': 2.20,
            't_half_life_d': 90,
            'ec50_ug_L': 100
        },
        'PFOS': {
            'M_product_%': 0.001,
            'M_conc_ug_L': 300,
            'log_koc': 3.65,
            't_half_life_d': 14965,
            'ec50_ug_L': 480
        },
        'Endosulfan': {
            'M_product_%': 0.005,
            'M_conc_ug_L': 20,
            'log_koc': 4.3,
            't_half_life_d': 60,
            'ec50_ug_L': 560
        }
    }

    # 2. Calculate aquifer properties
    # Based on mass balance in 1 m^3 of aquifer material
    # m_total = m_solid + m_water
    # V_total = V_solid + V_water = 1 m^3
    # m_water = 0.3 * m_total; m_solid = 0.7 * m_total
    # V_solid = m_solid / rho_solid; V_water = m_water / rho_water
    # 1 = (0.7*m_total/1346) + (0.3*m_total/1000)
    m_total_per_m3 = 1 / (0.7 / coal_density_kg_m3 + 0.3 / 1000)
    m_solid_per_m3 = 0.7 * m_total_per_m3 # This is the bulk density
    v_water_per_m3 = (0.3 * m_total_per_m3) / 1000 # This is the porosity

    rho_b = m_solid_per_m3
    n = v_water_per_m3
    
    # Calculate cross-sectional area and dispersion coefficient
    A = fracture_length_m * fracture_height_m  # Area perpendicular to flow
    D = dispersivity_m * advection_velocity_m_d # Dispersion coefficient m^2/d

    print("--- Aquifer and System Properties ---")
    print(f"Calculated Bulk Density (ρb): {rho_b:.2f} kg/m³")
    print(f"Calculated Porosity (n): {n:.3f}")
    print(f"Flow Cross-Sectional Area (A): {A:.0f} m²")
    print(f"Dispersion Coefficient (D): {D:.2f} m²/d")
    print("-" * 35 + "\n")

    # Function to calculate concentration using the standard 1D Advection-Dispersion equation
    def calculate_concentration(t, M_ug, R, k_d, x_m, A_m2, n_por, D_m2d, v_md):
        if t <= 0:
            return 0.0
        
        # Retarded velocity and dispersion
        v_retarded = v_md / R
        D_retarded = D_m2d / R

        # Denominator of main term
        denominator = A_m2 * n_por * math.sqrt(4 * math.pi * D_retarded * t)
        if denominator == 0: return 0.0

        # Numerator of exponential term
        exp_num = (x_m - v_retarded * t)**2
        exp_den = 4 * D_retarded * t
        if exp_den == 0: return 0.0

        try:
            conc_ug_m3 = (M_ug / denominator) * math.exp(-exp_num / exp_den) * math.exp(-k_d * t)
        except (ValueError, OverflowError):
            return 0.0
            
        return conc_ug_m3

    results = {}
    total_rq = 0

    for chemical, data in chemicals_data.items():
        print(f"--- Calculating for {chemical} ---")
        
        # Calculate initial mass
        M_ug = data['M_product_%'] * total_water_L * data['M_conc_ug_L']
        
        # Calculate chemical-specific parameters
        k = math.log(2) / data['t_half_life_d']
        koc = 10**data['log_koc']
        kd_L_kg = koc * coal_foc
        kd_m3_kg = kd_L_kg * 0.001 # Convert L to m^3
        R = 1 + (rho_b * kd_m3_kg) / n

        # Define the objective function for optimization (negative of concentration)
        objective_func = lambda t: -calculate_concentration(t, M_ug, R, k, distance_to_spring_m, A, n, D, advection_velocity_m_d)

        # Estimate peak arrival time to set optimization bounds
        t_peak_est = (distance_to_spring_m * R) / advection_velocity_m_d
        
        # Use SciPy to find the time of peak concentration
        # Bounds are set widely around the estimate to ensure the peak is found
        res = optimize.minimize_scalar(objective_func, bounds=(t_peak_est * 0.1, t_peak_est * 5.0), method='bounded')
        
        peak_time_days = res.x
        max_conc_ug_m3 = -res.fun
        max_conc_ug_L = max_conc_ug_m3 / 1000.0 # Convert from µg/m³ to µg/L

        rq = max_conc_ug_L / data['ec50_ug_L']
        total_rq += rq
        
        results[chemical] = {
            'peak_conc_ug_L': max_conc_ug_L,
            'peak_time_years': peak_time_days / 365.25,
            'RQ': rq
        }

        # Print the values used in the calculation for this chemical
        print("Parameters for final equation:")
        print(f"  Total Mass (M): {M_ug:,.0f} µg")
        print(f"  Retardation Factor (R): {R:,.2f}")
        print(f"  Decay Constant (k): {k:.6f} per day")
        print(f"  Peak Arrival Time (t): {peak_time_days:,.0f} days ({peak_time_days/365.25:,.1f} years)")
        print(f"Resulting Peak Concentration (C_max): {max_conc_ug_L:.6f} µg/L")
        print(f"Risk Quotient (RQ = C_max / EC50): {rq:.6f}")
        print("-" * 35 + "\n")

    # Determine the highest concentration
    highest_conc_chemical = max(results, key=lambda k: results[k]['peak_conc_ug_L'])
    highest_conc_value = results[highest_conc_chemical]['peak_conc_ug_L']

    # Final Summary
    print("--- Final Summary ---")
    print(f"The highest individual chemical concentration reaching the spring is for {highest_conc_chemical}.")
    print(f"Highest Concentration: {highest_conc_value:.6f} µg/L")
    print("\nAssessing the combined effect of the mixture using the Hazard Index (additive effects):")
    print(f"Hazard Index (HI = Σ RQ): {total_rq:.6f}")
    if total_rq >= 1:
        print("Effect: The Hazard Index is >= 1, indicating a potential risk from the additive effect of the chemical mixture.")
    else:
        print("Effect: The Hazard Index is < 1, suggesting the mixture is unlikely to pose a risk to the algal community.")
        
    return highest_conc_chemical, highest_conc_value, total_rq

if __name__ == '__main__':
    chem, conc, hi = solve_transport_problem()
    effect = "additive, with a potential risk" if hi >= 1 else "additive, with low potential for risk"
    print(f"<<<The highest concentration is {conc:.6f} μg/L for {chem}. The effect of the mixture is {effect}, with a Hazard Index of {hi:.6f}.>>>")
