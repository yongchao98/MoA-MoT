import math

def solve_fish_accumulation():
    """
    Solves the fish chemical accumulation problem based on the provided parameters.
    """
    # Step 1: Define constants and assumptions
    
    # Environment parameters
    V = 10000  # L (Volume of freshwater body)
    Qin = 900  # L/d (Inflow rate)
    # Assumption: To prevent the environment from drying out, assume Qout = Qin
    Qout = 900 # L/d (Outflow rate, was 1600)
    foc = 0.001  # 0.1% organic carbon
    t = 365  # days
    C0 = 0 # ng/L (Initial concentration in pristine water)

    # Fish parameters
    Mfish = 1000  # g (average weight)
    IRfood = 20  # g/day (ingestion rate)
    AFgills = 0.8  # absorption fraction through gills
    AFfood = 0.9  # absorption fraction for food
    Qgills = 100  # L/day (gill flow rate)

    # PFOS parameters
    pfos_params = {
        'name': 'PFOS',
        'Cin': 2.6,  # ng/L
        'half_life_years': 91,
        'logKow': 4.0,
        'Cfood': 100,  # ng/g
        'kelim': 0.069,  # days^-1
        'Cfish': 10  # ng/g
    }

    # PFOA parameters
    pfoa_params = {
        'name': 'PFOA',
        'Cin': 211300,  # ng/L
        'half_life_years': 238,
        'logKow': 4.5,
        'Cfood': 100,  # ng/g
        'kelim': 0.023,  # days^-1
        'Cfish': 10  # ng/g
    }
    
    # Calculate hydraulic flushing rate 'r'
    r = Qout / V  # days^-1

    print(f"Assumption: To ensure the water body does not dry up, the outflow rate (Qout) is set equal to the inflow rate (Qin) = {Qin} L/d.\n")

    total_accumulation = 0
    
    # Loop through each chemical to perform calculations
    for params in [pfos_params, pfoa_params]:
        print(f"--- Calculating for {params['name']} ---")
        
        # Step 2: Calculate Water Concentration C(t)
        # a) Calculate Koc
        logKoc = 0.81 * params['logKow'] + 0.01
        Koc = 10**logKoc
        
        # b) Calculate degradation rate constant k
        half_life_days = params['half_life_years'] * 365.25
        k = math.log(2) / half_life_days

        # c) Use Kd = Koc in the provided formula's denominator structure
        # Denominator of C_ss term
        denominator_ct = Qout + Qin * (1 + Koc * foc)

        # d) Calculate C(t=365)
        C_ss_term = (params['Cin'] * Qin) / denominator_ct
        exponent_term = math.exp(-(k + r) * t)
        C_t_365 = C_ss_term * (1 - exponent_term) + C0 * exponent_term

        print(f"Water concentration C({t} days) for {params['name']}: {C_t_365:.4f} ng/L")

        # Step 3: Calculate POP Accumulation Rate
        uptake_gills = C_t_365 * Qgills * AFgills
        uptake_food = params['Cfood'] * IRfood * AFfood
        elimination = params['kelim'] * params['Cfish'] * Mfish
        
        accumulation_rate = uptake_gills + uptake_food - elimination
        total_accumulation += accumulation_rate

        print("Equation: Accumulation = C(t) * Qgills * AFgills + Cfood * IRfood * AFfood - kelim * Cfish * Mfish")
        print(f"Accumulation Rate for {params['name']} = {C_t_365:.4f} * {Qgills} * {AFgills} + {params['Cfood']} * {IRfood} * {AFfood} - {params['kelim']} * {params['Cfish']} * {Mfish}")
        print(f"Accumulation Rate for {params['name']} = {accumulation_rate:.4f} ng/day\n")

    # Step 4: Final Total Accumulation
    print("--- Total Accumulation Rate ---")
    print(f"The total chemical accumulation rate is the sum of the rates for PFOS and PFOA.")
    print(f"Total Accumulation Rate = Accumulation(PFOS) + Accumulation(PFOA)")
    print(f"Total Accumulation Rate = {total_accumulation:.2f} ng/day")
    
    print(f"\n<<<Total chemical accumulation rate in fish is {total_accumulation:.2f} ng/day>>>")

solve_fish_accumulation()