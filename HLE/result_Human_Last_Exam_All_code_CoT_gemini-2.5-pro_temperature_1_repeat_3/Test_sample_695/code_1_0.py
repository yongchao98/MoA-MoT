import math

def calculate_pop_accumulation():
    """
    Calculates the net accumulation rate of PFOS and PFOA in fish
    based on the parameters provided in the problem description.
    """

    # --- Environmental and Fish Parameters (Shared) ---
    V_water = 10000  # L (Volume of freshwater body)
    Qin = 900      # L/d (Inflow rate)
    Qout = 1600    # L/d (Outflow rate)
    foc = 0.001    # dimensionless (0.1% organic carbon)
    t_days = 365   # days

    # Fish parameters
    M_fish = 1000  # g (Fish weight)
    C_food = 100   # ng/g (Concentration in food)
    IR_food = 20   # g/day (Ingestion rate of food)
    AF_gills = 0.8 # dimensionless (Absorption fraction, gills)
    AF_food = 0.9  # dimensionless (Absorption fraction, food)
    Q_gills = 100  # L/day (Gill flow rate)
    C_fish = 10    # ng/g (Concentration in fish at time t)

    # --- Chemical-Specific Parameters ---
    # PFOS
    Cin_pfos = 2.6       # ng/L
    t_half_pfos = 91     # years
    logKow_pfos = 4.0
    kelim_pfos = 0.069   # day⁻¹

    # PFOA
    Cin_pfoa = 211300    # ng/L
    t_half_pfoa = 238    # years
    logKow_pfoa = 4.5
    kelim_pfoa = 0.023   # day⁻¹

    # --- Calculation Function ---
    def calculate_for_chemical(name, Cin, t_half_years, logKow, kelim):
        print(f"--- Calculating for {name} ---")

        # 1. Calculate water concentration C(t) at t = 365 days
        t_half_days = t_half_years * 365
        k_decay = math.log(2) / t_half_days
        r_flush = Qout / V_water
        
        # Using formula: log Koc = 0.81 * log Kow + 0.01
        logKoc = 0.81 * logKow + 0.01
        Koc = 10**logKoc
        
        # Using formula: Kd = Koc * foc
        Kd = Koc * foc
        
        # Using water concentration formula provided:
        # C(t) = (Cin * Qin / (Qout + Qin * (1 + Kd * foc))) * (1 - e^-(k + r)t) + C0 * e^-(k + r)t
        # Assuming C0 (initial concentration) = 0 for a pristine environment.
        C0 = 0
        
        denominator = Qout + Qin * (1 + Kd * foc)
        rate_constant = k_decay + r_flush
        
        C_t = (Cin * Qin / denominator) * (1 - math.exp(-rate_constant * t_days)) + C0 * math.exp(-rate_constant * t_days)
        
        print(f"Water concentration C({t_days}) for {name}: {C_t:.4f} ng/L")

        # 2. Calculate net accumulation rate in fish at t = 365 days
        # Using formula: POP accumulation = C(t) * Qgills * AFgills + Cfood * IRfood * AFfood - kelim * Cfish * Mfish
        
        uptake_gills = C_t * Q_gills * AF_gills
        uptake_food = C_food * IR_food * AF_food
        elimination = kelim * C_fish * M_fish
        
        net_accumulation = uptake_gills + uptake_food - elimination
        
        print(f"\nFinal Equation for {name}:")
        print(f"Net Accumulation = C(t) * Qgills * AFgills + Cfood * IRfood * AFfood - kelim * Cfish * Mfish")
        print("Net Accumulation (ng/day) = "
              f"{C_t:.4f} * {Q_gills} * {AF_gills} + "
              f"{C_food} * {IR_food} * {AF_food} - "
              f"{kelim} * {C_fish} * {M_fish}")
        
        print(f"Net Accumulation (ng/day) = "
              f"{uptake_gills:.2f} (from gills) + "
              f"{uptake_food:.2f} (from food) - "
              f"{elimination:.2f} (elimination)")
        
        print(f"Net Accumulation Rate for {name}: {net_accumulation:.2f} ng/day\n")
        return net_accumulation

    # --- Execute Calculations ---
    accumulation_pfos = calculate_for_chemical("PFOS", Cin_pfos, t_half_pfos, logKow_pfos, kelim_pfos)
    accumulation_pfoa = calculate_for_chemical("PFOA", Cin_pfoa, t_half_pfoa, logKow_pfoa, kelim_pfoa)

    total_accumulation = accumulation_pfos + accumulation_pfoa
    print("--------------------------------------------------")
    print(f"Total Net Accumulation Rate (PFOS + PFOA): {total_accumulation:.2f} ng/day")
    print(f"<<<{total_accumulation:.2f}>>>")

# Run the calculation
calculate_pop_accumulation()