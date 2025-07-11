import math

def solve_accumulation():
    """
    Calculates the accumulation of PFOS and PFOA in fish after 365 days
    based on the provided environmental and toxicological parameters.
    """

    # --- Environmental Parameters ---
    V = 10000      # L (Volume of freshwater body)
    Qin = 900      # L/d (Inflow rate)
    Qout = 1600    # L/d (Outflow rate)
    foc = 0.001    # 0.1% organic carbon
    t = 365        # days

    # --- Fish Parameters (Shared) ---
    Mfish = 1000   # g (Fish weight)
    Cfood = 100    # ng/g (Food concentration)
    IRfood = 20    # g/day (Ingestion rate)
    AFgills = 0.8  # (Absorption fraction, gills)
    AFfood = 0.9   # (Absorption fraction, food)
    Qgills = 100   # L/day (Gill flow rate)
    Cfish_0 = 10   # ng/g (Initial fish concentration)

    # --- Chemical-Specific Parameters Dictionary ---
    chemicals = {
        "PFOS": {
            "Cin": 2.6,
            "logKow": 4.0,
            "kelim": 0.069
        },
        "PFOA": {
            "Cin": 211300,
            "logKow": 4.5,
            "kelim": 0.023
        }
    }
    
    results = {}

    for name, params in chemicals.items():
        print(f"--- Calculating for {name} ---")

        # Part 1: Calculate Water Concentration at day 365, C(365)
        # We assume water concentration reaches steady state (Css) as t is large.
        logKoc = 0.81 * params["logKow"] + 0.01
        Koc = 10**logKoc
        
        # We assume Kd in the provided formula refers to the calculated Koc value.
        Css_numerator = params["Cin"] * Qin
        Css_denominator = Qout + Qin * (1 + Koc * foc)
        C_water_365 = Css_numerator / Css_denominator
        
        # Part 2: Calculate Fish Concentration at day 365, C_fish(365)
        # Total uptake rate is from water and food.
        uptake_from_water = C_water_365 * Qgills * AFgills
        uptake_from_food = Cfood * IRfood * AFfood
        total_uptake_rate = uptake_from_water + uptake_from_food
        
        # Calculate steady-state fish concentration (C_fish_ss)
        C_fish_ss = total_uptake_rate / (params["kelim"] * Mfish)
        
        # Calculate transient fish concentration at t=365
        exp_term = math.exp(-params["kelim"] * t)
        C_fish_365 = C_fish_ss * (1 - exp_term) + Cfish_0 * exp_term

        # Part 3: Calculate the Net POP Accumulation Rate at day 365
        # Rate = Uptake_Water + Uptake_Food - Elimination_at_t=365
        elimination_rate_365 = params["kelim"] * C_fish_365 * Mfish
        net_accumulation_rate_365 = total_uptake_rate - elimination_rate_365

        # Print the breakdown of the final equation as requested
        print("Breakdown of the POP accumulation equation at day 365 (values are in ng/day):")
        final_eq = f"POP Accumulation Rate = (C(t) * Qgills * AFgills) + (Cfood * IRfood * AFfood) - (kelim * Cfish * Mfish)"
        print(final_eq)
        
        print(f"Result = ({C_water_365:.4f} * {Qgills} * {AFgills}) + ({Cfood} * {IRfood} * {AFfood}) - ({params['kelim']} * {C_fish_365:.4f} * {Mfish})")
        print(f"Result = {uptake_from_water:.2f} + {uptake_from_food:.2f} - {elimination_rate_365:.2f}")
        print(f"Net Accumulation Rate for {name} at day 365 = {net_accumulation_rate_365:.4f} ng/day\n")
        
        results[name] = net_accumulation_rate_365

    return results

if __name__ == '__main__':
    final_rates = solve_accumulation()
    # The final answer is the net accumulation rate calculated for each chemical.
    print(f"<<<PFOS Net Accumulation Rate: {final_rates['PFOS']:.4f} ng/day, PFOA Net Accumulation Rate: {final_rates['PFOA']:.4f} ng/day>>>")
