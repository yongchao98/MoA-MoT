import math

def solve_accumulation():
    """
    Calculates the accumulation rate of PFOS and PFOA in fish based on the provided environmental and biological data.
    """
    # --- Environment Parameters ---
    V = 10000  # Volume of freshwater body (L)
    Qin = 900  # Inflow rate (L/d)
    Qout = 1600 # Outflow rate (L/d)
    foc = 0.001 # Fraction of organic carbon (0.1%)
    C0 = 0      # Initial concentration in water (ng/L)
    t = 365     # Time (days)

    # --- Fish Parameters ---
    Mfish = 1000   # Mass of fish (g)
    Cfood = 100    # Concentration in food (ng/g)
    IRfood = 20    # Ingestion rate (g/day)
    AFgills = 0.8  # Absorption fraction for gills
    AFfood = 0.9   # Absorption fraction for food
    Qgills = 100   # Gill flow rate (L/day)
    Cfish = 10     # Initial/constant concentration in fish (ng/g)

    # --- Chemical Parameters ---
    chemicals = {
        "PFOS": {
            "Cin": 2.6,         # Inflow concentration (ng/L)
            "half_life": 91,    # years
            "log_Kow": 4.0,
            "kelim": 0.069      # Elimination rate constant (days⁻¹)
        },
        "PFOA": {
            "Cin": 211300,      # Inflow concentration (ng/L)
            "half_life": 238,   # years
            "log_Kow": 4.5,
            "kelim": 0.023      # Elimination rate constant (days⁻¹)
        }
    }

    # --- Calculations for each chemical ---
    results = {}
    for name, params in chemicals.items():
        # Step 1: Calculate water concentration C(t)
        
        # Degradation rate constant (k)
        half_life_days = params["half_life"] * 365.25
        k = math.log(2) / half_life_days

        # Organic carbon-water partition coefficient (Koc)
        log_Koc = 0.81 * params["log_Kow"] + 0.01
        Koc = 10**log_Koc
        
        # We will assume Kd (partition coefficient) is represented by Koc
        Kd = Koc

        # Hydraulic flushing rate (r)
        r = Qout / V

        # Loss rate sum
        loss_rate = k + r

        # Calculate steady-state concentration part of the equation
        numerator = params["Cin"] * Qin
        denominator = Qout + Qin * (1 + Kd * foc)
        Css_part = numerator / denominator

        # Calculate C(t) using the provided equation
        # C(t) = (Css_part) * (1 - e^-(k+r)t) + C0 * e^-(k+r)t
        # Since C0 is 0, the second term is zero.
        Ct = Css_part * (1 - math.exp(-loss_rate * t))
        
        # Step 2: Calculate POP accumulation rate in fish
        
        # Uptake from gills
        uptake_gills = Ct * Qgills * AFgills
        
        # Uptake from food
        uptake_food = Cfood * IRfood * AFfood
        
        # Elimination
        elimination = params["kelim"] * Cfish * Mfish
        
        # Net accumulation rate (ng/day)
        accumulation_rate = uptake_gills + uptake_food - elimination
        
        results[name] = {
            "Ct": Ct,
            "uptake_gills": uptake_gills,
            "uptake_food": uptake_food,
            "elimination": elimination,
            "accumulation_rate": accumulation_rate
        }

    # --- Print Results ---
    pfos_res = results["PFOS"]
    pfoa_res = results["PFOA"]

    print("--- PFOS Accumulation Rate Calculation ---")
    print(f"The concentration of PFOS in water after {t} days, C({t}), is {pfos_res['Ct']:.4f} ng/L.")
    print("POP accumulation = (C(t) * Qgills * AFgills) + (Cfood * IRfood * AFfood) - (kelim * Cfish * Mfish)")
    print(f"PFOS accumulation = ({pfos_res['Ct']:.4f} * {Qgills} * {AFgills}) + ({Cfood} * {IRfood} * {AFfood}) - ({chemicals['PFOS']['kelim']} * {Cfish} * {Mfish})")
    print(f"PFOS accumulation = {pfos_res['uptake_gills']:.2f} ng/day + {pfos_res['uptake_food']:.2f} ng/day - {pfos_res['elimination']:.2f} ng/day")
    print(f"Final PFOS accumulation rate = {pfos_res['accumulation_rate']:.2f} ng/day\n")

    print("--- PFOA Accumulation Rate Calculation ---")
    print(f"The concentration of PFOA in water after {t} days, C({t}), is {pfoa_res['Ct']:.4f} ng/L.")
    print("POP accumulation = (C(t) * Qgills * AFgills) + (Cfood * IRfood * AFfood) - (kelim * Cfish * Mfish)")
    print(f"PFOA accumulation = ({pfoa_res['Ct']:.4f} * {Qgills} * {AFgills}) + ({Cfood} * {IRfood} * {AFfood}) - ({chemicals['PFOA']['kelim']} * {Cfish} * {Mfish})")
    print(f"PFOA accumulation = {pfoa_res['uptake_gills']:.2f} ng/day + {pfoa_res['uptake_food']:.2f} ng/day - {pfoa_res['elimination']:.2f} ng/day")
    print(f"Final PFOA accumulation rate = {pfoa_res['accumulation_rate']:.2f} ng/day")
    
    # Final answer in the specified format
    final_answer = f"PFOS: {pfos_res['accumulation_rate']:.2f} ng/day, PFOA: {pfoa_res['accumulation_rate']:.2f} ng/day"
    print(f"\n<<<{final_answer}>>>")

solve_accumulation()