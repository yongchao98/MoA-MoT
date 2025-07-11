import math

def solve_accumulation():
    """
    Calculates the total chemical accumulation rate in fish for PFOS and PFOA
    based on the provided environmental and physiological parameters and formulas.
    """
    # --- Environmental and Fish Parameters ---
    V = 10000.0   # Volume of freshwater body (L)
    Qin = 900.0     # Inflow rate (L/d)
    Qout = 1600.0   # Outflow rate (L/d)
    foc = 0.001     # Organic carbon fraction (0.1%)
    t = 365.0       # Time (days)

    Mfish = 1000.0  # Weight of fish (g)
    Cfood = 100.0   # Food concentration of each chemical (ng/g)
    IRfood = 20.0   # Ingestion rate of food (g/day)
    AFgills = 0.8   # Gill absorption fraction
    AFfood = 0.9    # Food absorption fraction
    Qgills = 100.0  # Gill flow rate (L/day)
    Cfish_initial = 10.0 # Initial concentration in fish (ng/g)
    C0 = 0.0 # Initial concentration in pristine water (ng/L)

    # --- Chemical-specific Parameters ---
    # PFOS
    params_pfos = {
        'name': 'PFOS',
        'Cin': 2.6,
        'half_life_years': 91.0,
        'log_Kow': 4.0,
        'kelim': 0.069
    }
    # PFOA
    params_pfoa = {
        'name': 'PFOA',
        'Cin': 211300.0,
        'half_life_years': 238.0,
        'log_Kow': 4.5,
        'kelim': 0.023
    }

    total_accumulation_rate = 0
    
    # Hydraulic flushing rate is constant for the system
    r = Qout / V
    
    for params in [params_pfos, params_pfoa]:
        name = params['name']
        Cin = params['Cin']
        half_life_years = params['half_life_years']
        log_Kow = params['log_Kow']
        kelim = params['kelim']

        print(f"--- Calculating for {name} ---")

        # Step 1: Calculate C_water(t=365)
        half_life_days = half_life_years * 365.25
        k = math.log(2) / half_life_days
        log_Koc = 0.81 * log_Kow + 0.01
        Koc = 10**log_Koc
        
        # Using the provided C(t) formula and interpreting 'Kd' as 'Koc'
        denominator = Qout + Qin * (1 + Koc * foc)
        exp_term = math.exp(-(k + r) * t)
        C_water_365 = (Cin * Qin / denominator) * (1 - exp_term) + C0 * exp_term

        # Step 2: Calculate POP accumulation rate
        uptake_gills = C_water_365 * Qgills * AFgills
        uptake_food = Cfood * IRfood * AFfood
        elimination = kelim * Cfish_initial * Mfish
        accumulation_rate = uptake_gills + uptake_food - elimination
        
        total_accumulation_rate += accumulation_rate
        
        print(f"The 'POP accumulation' is interpreted as the instantaneous accumulation rate at t=365 days.")
        print(f"{name} Accumulation Rate = C(t) * Qgills * AFgills + Cfood * IRfood * AFfood - kelim * Cfish * Mfish")
        print(f"{name} Accumulation Rate = {C_water_365:.4f} ng/L * {Qgills} L/day * {AFgills} + {Cfood} ng/g * {IRfood} g/day * {AFfood} - {kelim} day^-1 * {Cfish_initial} ng/g * {Mfish} g")
        print(f"{name} Accumulation Rate = {accumulation_rate:.2f} ng/day\n")

    print("--- Total Accumulation Rate ---")
    print(f"The total chemical accumulation rate in the fish at day 365 is the sum of the rates for each chemical.")
    print(f"Total Accumulation Rate = {total_accumulation_rate:.2f} ng/day")
    
    # Final answer output
    print(f"<<<{total_accumulation_rate:.2f}>>>")

solve_accumulation()