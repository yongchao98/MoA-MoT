import math

def calculate_accumulation():
    """
    Calculates the accumulation rate of PFOS and PFOA in fish based on the provided environmental and chemical data.
    """
    # --- Given Parameters ---

    # Environment
    V_water = 10000  # L
    Qin = 900      # L/d
    Qout = 1600    # L/d
    foc = 0.001    # 0.1% organic carbon as a fraction

    # Time
    t = 365  # days

    # Fish properties
    M_fish = 1000   # g
    IR_food = 20    # g/day
    AF_gills = 0.8  # absorption fraction gills
    AF_food = 0.9   # absorption fraction food
    Q_gills = 100   # L/day
    C_food = 100    # ng/g of each chemical
    C_fish = 10     # ng/g of each chemical

    # Chemical properties dictionary
    chemicals = {
        "PFOS": {
            "Cin": 2.6,
            "half_life_years": 91,
            "log_Kow": 4.0,
            "kelim": 0.069
        },
        "PFOA": {
            "Cin": 211300,
            "half_life_years": 238,
            "log_Kow": 4.5,
            "kelim": 0.023
        }
    }

    print("Step 1: Calculating the chemical concentration in water, C(t), at t=365 days.\n")
    
    # --- Step 1: Calculate Water Concentration C(t) ---
    r = Qout / V_water  # Flushing rate in days^-1
    C0 = 0  # Assuming initial concentration in water is zero

    water_concentrations = {}
    
    for name, props in chemicals.items():
        # Convert half-life from years to days
        t_half_days = props["half_life_years"] * 365
        # Calculate degradation rate constant k
        k = math.log(2) / t_half_days

        # Calculate log Koc and then Koc
        log_Koc = 0.81 * props["log_Kow"] + 0.01
        Koc = 10**log_Koc

        # Calculate Kd = Koc * foc
        Kd = Koc * foc

        # Calculate steady-state concentration C_ss using the provided formula
        C_ss_numerator = props["Cin"] * Qin
        C_ss_denominator = Qout + Qin * (1 + Kd * foc)
        C_ss = C_ss_numerator / C_ss_denominator

        # Calculate the overall rate constant for the exponent
        rate_lambda = k + r

        # Calculate C(t) using the expression for concentration over time
        C_t = C_ss * (1 - math.exp(-rate_lambda * t)) + C0 * math.exp(-rate_lambda * t)
        
        water_concentrations[name] = C_t

    print("Step 2: Calculating the net accumulation rate in fish at day 365.\n")
    print("This is calculated using the formula:")
    print("POP accumulation = C(t) * Qgills * AFgills + Cfood * IRfood * AFfood - kelim * Cfish * Mfish\n")

    # --- Step 2: Calculate POP Accumulation Rate ---
    for name, props in chemicals.items():
        C_t_val = water_concentrations[name]
        kelim_val = props["kelim"]

        # Calculate each term of the accumulation equation
        uptake_gills = C_t_val * Q_gills * AF_gills
        uptake_food = C_food * IR_food * AF_food
        elimination = kelim_val * C_fish * M_fish

        # Calculate the final net accumulation rate
        net_accumulation_rate = uptake_gills + uptake_food - elimination

        # Print the final answer, showing the equation with all numbers plugged in
        print(f"--- For {name} ---")
        equation_str = (
            f"Net Accumulation Rate = "
            f"({C_t_val:.4f} ng/L * {Q_gills} L/day * {AF_gills}) "
            f"+ ({C_food} ng/g * {IR_food} g/day * {AF_food}) "
            f"- ({kelim_val} days⁻¹ * {C_fish} ng/g * {M_fish} g)"
        )
        print(equation_str)
        result_str = (
            f"                      = {uptake_gills:.2f} ng/day "
            f"+ {uptake_food:.2f} ng/day "
            f"- {elimination:.2f} ng/day "
            f"= {net_accumulation_rate:.2f} ng/day\n"
        )
        print(result_str)

calculate_accumulation()