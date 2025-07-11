import math

def solve_environmental_transport():
    """
    Calculates the peak concentration of three chemicals at a spring and assesses the mixture effect,
    based on a provided advection-diffusion equation.
    """
    # --- 1. Define initial parameters and constants ---

    # Environmental parameters from the problem description
    V_water_total = 1_500_000  # L
    coal_density = 1346  # kg/m^3 (Assumed as bulk density rho_b)
    foc = 0.50  # fraction of organic carbon
    theta = 0.30  # water content (porosity)
    distance_to_spring = 100  # m (x)
    advection_velocity = 0.5  # m/d (v)
    dispersivity = 0.5  # m (alpha_L)
    fracture_height = 10 # m
    fracture_width = 0.005 # m

    # Calculate parameters used in the equation for all chemicals
    # Dispersion coefficient D = alpha_L * v
    D = dispersivity * advection_velocity
    # Cross-sectional area A = height * width
    A = fracture_height * fracture_width

    # Define chemical-specific properties
    chemicals = {
        "Atrazine": {
            "product_percentage": 0.01, "concentration_in_product": 40, "log_koc": 2.20,
            "half_life_days": 90, "ec50": 100
        },
        "PFOS": {
            "product_percentage": 0.001, "concentration_in_product": 300, "log_koc": 3.65,
            "half_life_days": 14965, "ec50": 480
        },
        "Endosulfan": {
            "product_percentage": 0.005, "concentration_in_product": 20, "log_koc": 4.3,
            "half_life_days": 60, "ec50": 560
        }
    }

    # --- 2. Perform calculations for each chemical ---

    results = {}
    highest_conc = -1
    highest_conc_chem = None
    
    print("--- Calculation Steps & Results ---")
    print("This script solves for the peak concentration of each chemical at the spring (x=100m) using the provided Advection-Diffusion Equation.")
    print("Key Assumptions:")
    print("1. The term '1/(1+KdCgw)' is interpreted as '1/R', where R is the retardation factor: R = 1 + (bulk_density / porosity) * Kd.")
    print("2. The equation's output is treated as mass per length and converted to volumetric concentration by dividing by (fracture_area * porosity).\n")


    for name, props in chemicals.items():
        print(f"--- Calculating for {name} ---")

        # a. Calculate Total Mass (MTotal) in ug
        vol_product = V_water_total * props["product_percentage"]
        MTotal = vol_product * props["concentration_in_product"]
        print(f"Total Injected Mass (MTotal): {MTotal:,.0f} μg")

        # b. Calculate decay rate constant (k) in 1/d
        k = math.log(2) / props["half_life_days"]
        print(f"Decay rate (k): {k:.6f} d⁻¹")
        
        # c. Calculate partitioning coefficient (Kd) in L/kg
        Koc = 10**props["log_koc"]
        Kd = Koc * foc
        print(f"Partitioning Coefficient (Kd): {Kd:.2f} L/kg")

        # d. Calculate Retardation Factor (R) based on assumption 1
        R = 1 + (coal_density / theta) * Kd
        print(f"Retardation Factor (R): {R:,.2f}")

        # e. Find time to peak concentration (t_peak) by solving the quadratic equation: at^2 + bt + c = 0
        a_quad = (advection_velocity**2 / (4 * D)) + k
        b_quad = 0.5
        c_quad = -(distance_to_spring**2 / (4 * D))
        
        # Solve for t using the quadratic formula (we need the positive root)
        t_peak = (-b_quad + math.sqrt(b_quad**2 - 4 * a_quad * c_quad)) / (2 * a_quad)
        print(f"Time to Peak Concentration (t_peak): {t_peak:.2f} days")

        # f. Calculate the peak concentration C_max using the provided formula
        # The result from the formula, C_max_len, is in units of μg/m
        C_max_len = (MTotal / math.sqrt(4 * math.pi * D * t_peak)) * \
                    math.exp(-((distance_to_spring - advection_velocity * t_peak)**2) / (4 * D * t_peak)) * \
                    math.exp(-k * t_peak) * \
                    (1 / R)
        
        # Convert to volumetric concentration ug/m^3 by dividing by (Area * porosity)
        C_max_vol_m3 = C_max_len / (A * theta)

        # Convert to ug/L (1 m^3 = 1000 L)
        C_max_ug_L = C_max_vol_m3 / 1000

        results[name] = {
            'C_max_ug_L': C_max_ug_L, 'EC50': props['ec50'], 'MTotal': MTotal, 'k': k, 'R': R,
            't_peak': t_peak, 'D': D, 'v': advection_velocity, 'x': distance_to_spring,
            'A': A, 'theta': theta
        }
        
        if C_max_ug_L > highest_conc:
            highest_conc = C_max_ug_L
            highest_conc_chem = name

        print(f"Peak Concentration at Spring: {C_max_ug_L:.4e} μg/L\n")

    # --- 3. Determine Highest Concentration and Mixture Effect ---

    print("--- Final Summary ---")
    print(f"The chemical with the highest concentration reaching the spring is {highest_conc_chem}.")

    # Display the full equation for the highest concentration chemical as requested
    chem_data = results[highest_conc_chem]
    print("\nFinal equation for the highest concentration chemical ({name}):".format(name=highest_conc_chem))
    equation_str = (
        f"C_max(μg/L) = ( {chem_data['MTotal']:.0f} / sqrt(4*π*{chem_data['D']}*{chem_data['t_peak']:.2f}) ) * "
        f"exp( -({chem_data['x']}-{chem_data['v']}*{chem_data['t_peak']:.2f})² / (4*{chem_data['D']}*{chem_data['t_peak']:.2f}) ) * "
        f"exp( -{chem_data['k']:.6f}*{chem_data['t_peak']:.2f} ) * (1 / {chem_data['R']:.1f}) / "
        f"({chem_data['A']}*{chem_data['theta']}) / 1000"
    )
    print(equation_str)
    print(f"Highest Concentration = {highest_conc:.4e} μg/L")

    # Mixture effect analysis
    print("\nMixture Effect Analysis:")
    RI = 0
    for name, res in results.items():
        RQ = res['C_max_ug_L'] / res['EC50']
        print(f"Risk Quotient (RQ) for {name}: {RQ:.2e}")
        RI += RQ

    print(f"\nTotal Risk Index (RI = Sum of RQs): {RI:.2e}")
    print("The Risk Index is calculated by summing the individual Risk Quotients. This 'Concentration Addition' approach assumes an ADDITIVE effect of the chemicals in the mixture.")

    final_answer_text = f"The highest concentration is {highest_conc:.4e} µg/L ({highest_conc_chem}), and the mixture effect is Additive."
    
    print(f"\n<<<{final_answer_text}>>>")

if __name__ == '__main__':
    solve_environmental_transport()