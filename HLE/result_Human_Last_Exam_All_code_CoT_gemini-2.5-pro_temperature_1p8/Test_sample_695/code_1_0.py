import math

def solve_bioaccumulation():
    # --- Environment Parameters ---
    V = 10000  # L (Volume of freshwater body)
    Qin = 900  # L/d (Inflow rate)
    Qout = 1600 # L/d (Outflow rate)
    foc = 0.001 # 0.1% organic carbon, as a fraction
    t = 365 # days

    # --- Fish Parameters ---
    Mfish = 1000  # g (Weight of fish)
    Cfood = 100  # ng/g (Concentration in food)
    IRfood = 20  # g/day (Ingestion rate)
    AFgills = 0.8  # Absorption fraction for gills
    AFfood = 0.9  # Absorption fraction for food
    Qgills = 100  # L/day (Gill flow rate)
    # Cfish_initial = 10 ng/g (not needed for steady-state calculation)

    # --- Chemical-Specific Parameters ---
    chemicals = {
        "PFOS": {
            "Cin": 2.6, # ng/L
            "t_half_years": 91,
            "logKow": 4.0,
            "kelim": 0.069 # days⁻¹
        },
        "PFOA": {
            "Cin": 211300, # ng/L
            "t_half_years": 238,
            "logKow": 4.5,
            "kelim": 0.023 # days⁻¹
        }
    }
    
    results = {}

    print("--- Calculation of Chemical Concentrations in Fish ---\n")

    for name, params in chemicals.items():
        print(f"--- For {name} ---")

        # Step 1: Calculate steady-state water concentration (C_water_ss)
        logKow = params["logKow"]
        t_half_days = params["t_half_years"] * 365
        Cin = params["Cin"]

        logKoc = 0.81 * logKow + 0.01
        Koc = 10**logKoc
        Kd = Koc * foc
        
        # Per the problem description, using the provided (dimensionally unusual) formula.
        # C_ss = (Cin * Qin) / (Qout + Qin * (1 + Kd * foc))
        numerator_c_water = Cin * Qin
        denominator_c_water = Qout + Qin * (1 + Kd * foc)
        C_water_ss = numerator_c_water / denominator_c_water

        print("1. Calculate steady-state water concentration (C_water_ss):")
        print(f"   logKoc = 0.81 * {logKow} + 0.01 = {logKoc:.4f}")
        print(f"   Koc = 10^{logKoc:.4f} = {Koc:.4f}")
        print(f"   Kd = {Koc:.4f} * {foc} = {Kd:.4f}")
        print(f"   C_water_ss = ({Cin} * {Qin}) / ({Qout} + {Qin} * (1 + {Kd:.4f} * {foc}))")
        print(f"   C_water_ss = {numerator_c_water} / {denominator_c_water:.4f} = {C_water_ss:.4f} ng/L\n")

        # Step 2: Calculate steady-state fish concentration (C_fish_ss)
        kelim = params["kelim"]

        # C_fish_ss = (Uptake_gills + Uptake_food) / Elimination_term
        # C_fish_ss = (C_water_ss * Qgills * AFgills + Cfood * IRfood * AFfood) / (kelim * Mfish)
        uptake_gills = C_water_ss * Qgills * AFgills
        uptake_food = Cfood * IRfood * AFfood
        elim_denominator = kelim * Mfish
        
        C_fish_ss = (uptake_gills + uptake_food) / elim_denominator
        results[name] = C_fish_ss

        print("2. Calculate steady-state fish concentration (C_fish_ss):")
        print("   Formula: C_fish_ss = (C_water_ss * Qgills * AFgills + Cfood * IRfood * AFfood) / (kelim * Mfish)")
        print(f"   C_fish_ss = (({C_water_ss:.4f} * {Qgills} * {AFgills}) + ({Cfood} * {IRfood} * {AFfood})) / ({kelim} * {Mfish})")
        print(f"   C_fish_ss = ({uptake_gills:.2f} + {uptake_food:.2f}) / {elim_denominator:.2f}")
        print(f"   C_fish_ss = {(uptake_gills + uptake_food):.2f} / {elim_denominator:.2f}")
        print(f"   Final C_fish_ss = {C_fish_ss:.4f} ng/g\n")

    pfos_result = results['PFOS']
    pfoa_result = results['PFOA']
    final_answer = f"PFOS concentration: {pfos_result:.2f} ng/g, PFOA concentration: {pfoa_result:.2f} ng/g"
    print(f"\n--- Final Answer ---\n{final_answer}")
    print(f"\n<<<PFOS concentration: {pfos_result:.2f} ng/g, PFOA concentration: {pfoa_result:.2f} ng/g>>>")

solve_bioaccumulation()