import math

def calculate_fish_concentration():
    """
    This function calculates the concentration of PFOS and PFOA in fish
    after 365 days based on the environmental parameters provided.
    """
    # --- Shared Environmental and Fish Parameters ---
    Qin = 900      # Inflow rate (L/d)
    Qout = 1600    # Outflow rate (L/d)
    foc = 0.001    # Fraction of organic carbon (0.1%)
    t = 365        # Time (days)

    Mfish = 1000   # Mass of fish (g)
    Cfood = 100    # Concentration in food (ng/g)
    IRfood = 20    # Ingestion rate of food (g/day)
    AFgills = 0.8  # Absorption fraction from gills
    AFfood = 0.9   # Absorption fraction from food
    Qgills = 100   # Gill flow rate (L/day)
    Cfish_initial = 10  # Initial concentration in fish (ng/g)

    # --- Chemical-Specific Parameters ---
    chemicals = {
        "PFOS": {
            "Cin": 2.6,         # Inflow concentration (ng/L)
            "logKow": 4.0,      # Log Octanol-Water Partition Coefficient
            "kelim": 0.069      # Elimination rate constant (days⁻¹)
        },
        "PFOA": {
            "Cin": 211300,      # Inflow concentration (ng/L)
            "logKow": 4.5,      # Log Octanol-Water Partition Coefficient
            "kelim": 0.023      # Elimination rate constant (days⁻¹)
        }
    }
    
    # Store results to present at the end
    results = {}

    print("--- Calculating Final Fish Concentrations after 365 days ---\n")

    for name, params in chemicals.items():
        # --- Step 1: Calculate Water Concentration C(t) ---
        # Using log Koc = 0.81 * log Kow + 0.01 and Kd = Koc * foc
        logKoc = 0.81 * params["logKow"] + 0.01
        Koc = 10**logKoc
        Kd = Koc * foc
        
        # Using the provided C(t) equation, which approximates to steady-state C_water
        # C_water = (Cin * Qin) / (Qout + Qin * (1 + Kd * foc))
        # The prompt formula has (1 + Kd * foc), which we will follow exactly.
        denominator = Qout + Qin * (1 + Kd * foc)
        C_water = (params["Cin"] * Qin) / denominator

        # --- Step 2: Determine Fish Uptake Parameters ---
        # Uptake Rate = Uptake from water + Uptake from food
        uptake_water = C_water * Qgills * AFgills
        uptake_food = Cfood * IRfood * AFfood
        
        # Total uptake rate constant K ((ng/g)/day)
        K = (uptake_water + uptake_food) / Mfish

        # --- Step 3 & 4: Solve for and Calculate Final Fish Concentration ---
        # The solution to the differential equation for fish concentration is:
        # C_fish(t) = (K / kelim) * (1 - e^(-kelim * t)) + C_fish(0) * e^(-kelim * t)
        kelim = params["kelim"]
        steady_state_fish_conc = K / kelim
        exp_term = math.exp(-kelim * t)
        
        final_fish_conc = steady_state_fish_conc * (1 - exp_term) + Cfish_initial * exp_term
        results[name] = final_fish_conc

        # --- Output the Final Equation with Numbers ---
        print(f"Final concentration equation for {name}:")
        print(f"C_fish({t}) = (K / k_elim) * (1 - e^(-k_elim * t)) + C_fish_initial * e^(-k_elim * t)")
        print(f"Calculated K = {K:.4f} (ng/g)/day")
        print(f"C_fish({t}) = ({K:.4f} / {kelim}) * (1 - e^(-{kelim} * {t})) + {Cfish_initial} * e^(-{kelim} * {t})")
        print(f"Resulting {name} concentration in fish = {final_fish_conc:.4f} ng/g\n")
        
    return results

if __name__ == '__main__':
    final_concentrations = calculate_fish_concentration()
    pfos_result = final_concentrations.get("PFOS", 0)
    pfoa_result = final_concentrations.get("PFOA", 0)
    # The final answer is wrapped in <<<>>>
    print(f"<<<PFOS: {pfos_result:.4f} ng/g, PFOA: {pfoa_result:.4f} ng/g>>>")
