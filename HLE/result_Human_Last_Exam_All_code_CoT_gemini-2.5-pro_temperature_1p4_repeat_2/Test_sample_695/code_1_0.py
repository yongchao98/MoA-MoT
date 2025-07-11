import math

def solve_fish_accumulation():
    """
    Calculates the net accumulation rate of PFOS and PFOA in fish after 365 days.
    """
    # --- Given Parameters ---

    # Environment
    V = 10000      # L (Volume of water body)
    Qin = 900      # L/d (Inflow rate)
    Qout = 1600    # L/d (Outflow rate)
    foc = 0.001    # 0.1% organic carbon, as a fraction
    C0 = 0         # ng/L (Initial water concentration, assumed pristine)
    t = 365        # days

    # Fish
    Mfish = 1000   # g (Fish weight)
    Cfood = 100    # ng/g (Concentration in food)
    IRfood = 20    # g/day (Ingestion rate)
    AFgills = 0.8  # Absorption fraction, gills
    AFfood = 0.9   # Absorption fraction, food
    Qgills = 100   # L/day (Gill flow rate)

    # Chemical-specific parameters
    chemicals = {
        "PFOS": {
            "Cin": 2.6,        # ng/L
            "t_half_years": 91,
            "logKow": 4.0,
            "kelim": 0.069,    # day⁻¹
            "Cfish": 10        # ng/g
        },
        "PFOA": {
            "Cin": 211300,     # ng/L
            "t_half_years": 238,
            "logKow": 4.5,
            "kelim": 0.023,    # day⁻¹
            "Cfish": 10        # ng/g
        }
    }

    # --- Calculations ---
    # Calculate hydraulic flushing rate (r), common for both
    r = Qout / V

    for name, params in chemicals.items():
        print(f"--- Calculations for {name} ---")

        # 1. Calculate Water Concentration C(t) at t=365 days
        
        # Convert half-life from years to days
        t_half_days = params["t_half_years"] * 365.25
        # Calculate degradation rate constant (k)
        k = math.log(2) / t_half_days
        # Calculate logKoc and Koc
        logKoc = 0.81 * params["logKow"] + 0.01
        Koc = 10**logKoc
        # Calculate Kd
        Kd = Koc * foc

        # Using the provided equation for C(t):
        denominator = Qout + Qin * (1 + Kd * foc)
        exponent_val = -(k + r) * t
        exp_term = math.exp(exponent_val)
        
        C_t = (params["Cin"] * Qin / denominator) * (1 - exp_term) + C0 * exp_term
        print(f"Water concentration of {name} at {t} days, C({t}): {C_t:.4f} ng/L")

        # 2. Calculate Net Accumulation Rate in fish at t=365 days
        uptake_gills = C_t * Qgills * AFgills
        uptake_food = Cfood * IRfood * AFfood
        elimination = params["kelim"] * params["Cfish"] * Mfish
        
        net_accumulation_rate = uptake_gills + uptake_food - elimination

        print(f"\nCalculating the Net Accumulation Rate for {name}:")
        print(f"Net Rate = (C_water * Qgills * AFgills) + (Cfood * IRfood * AFfood) - (kelim * Cfish * Mfish)")
        print(f"Net Rate = ({C_t:.4f} * {Qgills} * {AFgills}) + ({Cfood} * {IRfood} * {AFfood}) - ({params['kelim']} * {params['Cfish']} * {Mfish})")
        print(f"Net Rate = ({uptake_gills:.2f}) + ({uptake_food:.2f}) - ({elimination:.2f})")
        print(f"Final Net Accumulation Rate for {name} at Day {t}: {net_accumulation_rate:.2f} ng/day\n")

solve_fish_accumulation()