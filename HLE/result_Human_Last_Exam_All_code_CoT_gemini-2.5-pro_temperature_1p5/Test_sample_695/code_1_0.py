import math

def calculate_fish_accumulation():
    """
    Calculates the accumulation rate of PFOS and PFOA in fish based on environmental and chemical properties.
    """
    # --- General Environmental and Fish Parameters ---
    V = 10000  # Volume of freshwater body (L)
    Qin = 900  # Inflow rate (L/d)
    Qout = 1600 # Outflow rate (L/d)
    foc = 0.001 # Fraction of organic carbon (0.1%)
    t = 365    # Time (days)
    C0 = 0     # Initial concentration in environment (ng/L), assumed 0 for pristine environment

    # Fish parameters
    Mfish = 1000   # Mass of fish (g)
    Cfood = 100    # Concentration of each chemical in food (ng/g)
    IRfood = 20    # Ingestion rate of food (g/day)
    AFfood = 0.9   # Absorption fraction for food
    Qgills = 100   # Gill flow rate (L/day)
    AFgills = 0.8  # Absorption fraction through gills
    Cfish = 10     # Given concentration of chemicals in fish (ng/g)

    # --- Chemical-Specific Parameters ---
    chemicals = {
        "PFOS": {
            "Cin": 2.6,         # Inflow concentration (ng/L)
            "half_life_years": 91,
            "logKow": 4.0,
            "kelim": 0.069      # Elimination rate constant (days⁻¹)
        },
        "PFOA": {
            "Cin": 211300,
            "half_life_years": 238,
            "logKow": 4.5,
            "kelim": 0.023
        }
    }

    total_accumulation_rate = 0

    print("--- Starting Calculation ---")

    # --- Step 1: Calculate shared environmental rates ---
    # Hydraulic residence rate (r)
    r = Qout / V

    for name, params in chemicals.items():
        print(f"\n--- Calculating for {name} ---")

        # --- Step 2: Calculate Water Concentration C(t) ---
        # Degradation rate constant (k) from half-life
        half_life_days = params["half_life_years"] * 365.25
        k = math.log(2) / half_life_days

        # Partition coefficients (Koc and Kd)
        logKoc = 0.81 * params["logKow"] + 0.01
        Koc = 10**logKoc
        Kd = Koc * foc # The problem asks to use this value for Kd

        # Denominator from the provided C(t) equation.
        # Note: The term `(1 + Kd * foc)` seems unusual as Kd is already a function of foc.
        # However, we will proceed as per the provided formula.
        denominator = Qout + Qin * (1 + Kd * foc)
        
        # Exponent term for C(t) equation
        exponent_val = -(k + r) * t
        
        # Calculate water concentration C(t) at t=365 days
        C_t = (params["Cin"] * Qin / denominator) * (1 - math.exp(exponent_val)) + C0 * math.exp(exponent_val)

        print(f"Water concentration C(t) at {t} days: {C_t:.4f} ng/L")
        
        # --- Step 3: Calculate Net Chemical Accumulation Rate in Fish ---
        # The formula is: Accumulation = (Gill Uptake) + (Food Uptake) - (Elimination)

        # Gill Uptake = C(t) * Qgills * AFgills
        gill_uptake = C_t * Qgills * AFgills
        
        # Food Uptake = Cfood * IRfood * AFfood
        food_uptake = Cfood * IRfood * AFfood

        # Elimination = kelim * Cfish * Mfish
        elimination_loss = params["kelim"] * Cfish * Mfish
        
        # Net accumulation rate for the chemical
        net_accumulation = gill_uptake + food_uptake - elimination_loss
        total_accumulation_rate += net_accumulation

        # Print the equation with all the numbers
        print("\nCalculating accumulation rate using the formula:")
        print("Accumulation = C(t) * Qgills * AFgills + Cfood * IRfood * AFfood - kelim * Cfish * Mfish")
        print("Substituted values:")
        print(f"Accumulation = {C_t:.4f} ng/L * {Qgills} L/day * {AFgills} + {Cfood} ng/g * {IRfood} g/day * {AFfood} - {params['kelim']} 1/day * {Cfish} ng/g * {Mfish} g")
        
        print("\nBreakdown of rates:")
        print(f"  Gill Uptake: {gill_uptake:.2f} ng/day")
        print(f"  Food Uptake: {food_uptake:.2f} ng/day")
        print(f"  Elimination: {elimination_loss:.2f} ng/day")
        
        print(f"Net {name} accumulation rate: {net_accumulation:.2f} ng/day")

    print("\n--- Final Result ---")
    print(f"Total chemical accumulation rate (PFOS + PFOA) in fish at day 365 is: {total_accumulation_rate:.2f} ng/day")
    return total_accumulation_rate

if __name__ == '__main__':
    final_answer = calculate_fish_accumulation()
    print(f"\n<<<Total accumulation rate: {final_answer:.2f}>>>")
