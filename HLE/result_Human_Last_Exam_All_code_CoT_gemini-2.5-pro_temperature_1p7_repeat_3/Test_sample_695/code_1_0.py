import math

def calculate_fish_accumulation():
    """
    Calculates the accumulated concentration of PFOS and PFOA in fish over 365 days.
    """
    # --- Step 1: Define all given parameters ---
    # Environment parameters
    V_water = 10000  # L (Volume of freshwater body)
    Qin = 900  # L/d (Inflow rate)
    Qout = 1600  # L/d (Outflow rate)
    foc = 0.001  # 0.1% organic carbon
    t = 365  # days
    C0_water = 0 # ng/L (Initial concentration in pristine environment)

    # Fish parameters
    M_fish = 1000  # g (Average weight of fish)
    C_food = 100  # ng/g (Concentration in food)
    IR_food = 20  # g/day (Ingestion rate)
    AF_gills = 0.8  # Absorption fraction through gills
    AF_food = 0.9  # Absorption fraction for food
    Q_gills = 100  # L/day (Gill flow rate)
    C_fish_initial = 10 # ng/g (Initial concentration in fish)

    # Chemical-specific parameters
    # PFOS
    Cin_pfos = 2.6  # ng/L
    half_life_pfos_years = 91
    log_Kow_pfos = 4.0
    kelim_pfos = 0.069  # days⁻¹
    
    # PFOA
    Cin_pfoa = 211300  # ng/L
    half_life_pfoa_years = 238
    log_Kow_pfoa = 4.5
    kelim_pfoa = 0.023  # days⁻¹

    # --- Calculations for each chemical ---
    chemicals = {
        "PFOS": {
            "Cin": Cin_pfos,
            "half_life_years": half_life_pfos_years,
            "log_Kow": log_Kow_pfos,
            "kelim": kelim_pfos
        },
        "PFOA": {
            "Cin": Cin_pfoa,
            "half_life_years": half_life_pfoa_years,
            "log_Kow": log_Kow_pfoa,
            "kelim": kelim_pfoa
        }
    }
    
    print("--- Calculation of Chemical Accumulation in Fish over 365 Days ---\n")

    for name, params in chemicals.items():
        # --- Step 2 & 3: Calculate final water concentration C(t) ---
        
        # Convert half-life from years to days
        half_life_days = params["half_life_years"] * 365.25
        # Calculate degradation rate constant k
        k = math.log(2) / half_life_days
        
        # Calculate log Koc and Koc
        log_Koc = 0.81 * params["log_Kow"] + 0.01
        Koc = 10**log_Koc
        
        # Calculate Kd
        Kd = Koc * foc
        
        # Calculate residence rate r
        r = Qout / V_water
        
        # Calculate water concentration C(t) for t=365 days
        # Denominator of the main term in C(t) equation
        # Based on the provided formula: Qout + Qin * (1 + Kd * foc)
        C_t_denominator = Qout + Qin * (1 + Kd * foc)
        # Exponent term
        exponent = (k + r) * t
        
        C_t_365 = (params["Cin"] * Qin / C_t_denominator) * (1 - math.exp(-exponent)) + C0_water * math.exp(-exponent)

        # --- Step 4 & 5: Calculate net accumulation in fish ---
        uptake_gills = C_t_365 * Q_gills * AF_gills
        uptake_food = C_food * IR_food * AF_food
        elimination = params["kelim"] * C_fish_initial * M_fish
        
        net_accumulation_rate = uptake_gills + uptake_food - elimination
        total_accumulated_mass = net_accumulation_rate * t
        
        # --- Step 6: Print the final results ---
        print(f"--- {name} Accumulation ---")
        print(f"The final water concentration of {name} after {t} days is: {C_t_365:.4f} ng/L")
        print("The total accumulated mass is calculated using the equation:")
        print("Accumulated Mass = (C_water * Qgills * AFgills + Cfood * IRfood * AFfood - kelim * Cfish * Mfish) * days\n")
        print("Substituting the values:")
        print(f"Accumulated {name} (ng) = \n"
              f"  ( {C_t_365:.4f} ng/L * {Q_gills} L/day * {AF_gills}  (from gills)\n"
              f"  + {C_food} ng/g * {IR_food} g/day * {AF_food}    (from food)\n"
              f"  - {params['kelim']} days⁻¹ * {C_fish_initial} ng/g * {M_fish} g )   (elimination)\n"
              f"  * {t} days\n")
        print(f"Daily Net Accumulation Rate for {name}:")
        print(f"  ({uptake_gills:.2f} + {uptake_food:.2f} - {elimination:.2f}) ng/day = {net_accumulation_rate:.2f} ng/day\n")
        print(f"Total Accumulated Mass of {name} in fish over {t} days: {total_accumulated_mass:,.2f} ng")
        print("-" * 50 + "\n")

if __name__ == '__main__':
    calculate_fish_accumulation()
    # To directly output the final answer for the system, let's re-run a simplified calculation.
    # We need to present both PFOS and PFOA results. The prompt asks for a single answer format.
    # I will output the final accumulated mass for both, separated by a comma.

    # --- PFOS ---
    half_life_pfos_days = 91 * 365.25
    k_pfos = math.log(2) / half_life_pfos_days
    log_Koc_pfos = 0.81 * 4.0 + 0.01
    Koc_pfos = 10**log_Koc_pfos
    Kd_pfos = Koc_pfos * 0.001
    r_pfos = 1600 / 10000
    C_t_denom_pfos = 1600 + 900 * (1 + Kd_pfos * 0.001)
    exp_pfos = (k_pfos + r_pfos) * 365
    C_t_365_pfos = (2.6 * 900 / C_t_denom_pfos) * (1 - math.exp(-exp_pfos))
    uptake_gills_pfos = C_t_365_pfos * 100 * 0.8
    uptake_food_pfos = 100 * 20 * 0.9
    elimination_pfos = 0.069 * 10 * 1000
    net_rate_pfos = uptake_gills_pfos + uptake_food_pfos - elimination_pfos
    total_accumulated_pfos = net_rate_pfos * 365
    
    # --- PFOA ---
    half_life_pfoa_days = 238 * 365.25
    k_pfoa = math.log(2) / half_life_pfoa_days
    log_Koc_pfoa = 0.81 * 4.5 + 0.01
    Koc_pfoa = 10**log_Koc_pfoa
    Kd_pfoa = Koc_pfoa * 0.001
    r_pfoa = 1600 / 10000
    C_t_denom_pfoa = 1600 + 900 * (1 + Kd_pfoa * 0.001)
    exp_pfoa = (k_pfoa + r_pfoa) * 365
    C_t_365_pfoa = (211300 * 900 / C_t_denom_pfoa) * (1 - math.exp(-exp_pfoa))
    uptake_gills_pfoa = C_t_365_pfoa * 100 * 0.8
    uptake_food_pfoa = 100 * 20 * 0.9
    elimination_pfoa = 0.023 * 10 * 1000
    net_rate_pfoa = uptake_gills_pfoa + uptake_food_pfoa - elimination_pfoa
    total_accumulated_pfoa = net_rate_pfoa * 365
    
    # Format the final answer string
    answer_string = f"PFOS: {total_accumulated_pfos:,.2f} ng, PFOA: {total_accumulated_pfoa:,.2f} ng"
    print(f"<<<{answer_string}>>>")