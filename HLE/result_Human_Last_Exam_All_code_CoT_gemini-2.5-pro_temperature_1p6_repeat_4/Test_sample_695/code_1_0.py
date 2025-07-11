import math

def calculate_accumulation():
    """
    Calculates the accumulation of PFOS and PFOA in fish based on the provided environmental parameters and equations.
    """
    # --- Environment and Fish Parameters (constant for both chemicals) ---
    V_water = 10000  # L (Volume of freshwater body)
    Qin = 900      # L/d (Inflow rate)
    Qout = 1600    # L/d (Outflow rate)
    foc = 0.001    # dimensionless (Organic carbon fraction, 0.1%)
    t = 365        # days (Time period)
    C0 = 0         # ng/L (Initial water concentration, assumed 0 for pristine environment)

    M_fish = 1000  # g (Weight of fish)
    C_food = 100   # ng/g (Concentration in food)
    IR_food = 20   # g/day (Ingestion rate of food)
    AF_gills = 0.8 # dimensionless (Absorption fraction for gills)
    AF_food = 0.9  # dimensionless (Absorption fraction for food)
    Q_gills = 100  # L/day (Gill flow rate)
    C_fish = 10    # ng/g (Initial/current concentration in fish for elimination calc)
    
    # --- Chemical-specific parameters ---
    chemicals = {
        "PFOS": {
            "Cin": 2.6,         # ng/L (Inflow concentration)
            "half_life_y": 91,  # years
            "log_Kow": 4.0,     # log unit
            "kelim": 0.069      # days⁻¹ (Elimination rate constant)
        },
        "PFOA": {
            "Cin": 211300,      # ng/L
            "half_life_y": 238, # years
            "log_Kow": 4.5,     # log unit
            "kelim": 0.023      # days⁻¹
        }
    }
    
    # Hydraulic residence rate (same for both)
    r = Qout / V_water

    results = {}

    for name, params in chemicals.items():
        print(f"--- Calculating for {name} ---")

        # Step 1: Calculate parameters for C(t)
        half_life_d = params["half_life_y"] * 365
        k = math.log(2) / half_life_d
        
        log_Koc = 0.81 * params["log_Kow"] + 0.01
        Koc = math.pow(10, log_Koc)
        # Assuming Kd = Koc * foc, as is standard practice
        Kd = Koc * foc
        
        # Step 2: Calculate water concentration C(t) using the provided formula
        # C(t) = (Cin * Qin / (Qout + Qin * (1 + Kd * foc))) * (1 - e^-(k + r)t) + C0 * e^-(k + r)t
        # Since C0 = 0, the second term is zero.
        rate_term = k + r
        C_ss_numerator = params["Cin"] * Qin
        C_ss_denominator = Qout + Qin * (1 + Kd * foc)
        C_ss = C_ss_numerator / C_ss_denominator
        
        C_t = C_ss * (1 - math.exp(-rate_term * t))

        # Step 3: Calculate the terms for the POP accumulation equation
        # POP accumulation = C(t) * Qgills * AFgills + Cfood * IRfood * AFfood - kelim * Cfish * Mfish
        uptake_gills = C_t * Q_gills * AF_gills
        uptake_food = C_food * IR_food * AF_food
        elimination = params["kelim"] * C_fish * M_fish

        # Step 4: Calculate final accumulation rate
        accumulation = uptake_gills + uptake_food - elimination
        results[name] = accumulation
        
        print("\nCalculating accumulation rate using the formula:")
        print("POP_accumulation = C(t) * Q_gills * AF_gills + C_food * IR_food * AF_food - k_elim * C_fish * M_fish")
        print("\nPlugging in the numbers for {}:".format(name))
        print(f"Accumulation = {C_t:.4f} * {Q_gills} * {AF_gills} + {C_food} * {IR_food} * {AF_food} - {params['kelim']} * {C_fish} * {M_fish}")
        print(f"Accumulation = {uptake_gills:.2f} + {uptake_food:.2f} - {elimination:.2f}")
        print(f"\nThe net accumulation rate for {name} after {t} days is: {accumulation:.2f} ng/day\n")

    pfos_result = results.get("PFOS", 0)
    pfoa_result = results.get("PFOA", 0)
    
    # Return final answer in the specified format
    print(f"<<<PFOS accumulation: {pfos_result:.2f} ng/day, PFOA accumulation: {pfoa_result:.2f} ng/day>>>")


if __name__ == "__main__":
    calculate_accumulation()