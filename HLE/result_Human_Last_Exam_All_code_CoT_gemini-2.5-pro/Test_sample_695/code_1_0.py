import math

def solve_pop_accumulation():
    """
    Calculates the total accumulated concentration of PFOS and PFOA in fish
    after 365 days based on the provided environmental and chemical parameters.
    """

    # --- Given Parameters ---

    # Environment
    V = 10000      # L (Volume)
    Qin = 900      # L/d (Inflow rate)
    Qout = 1600    # L/d (Outflow rate)
    foc = 0.001    # 0.1% (fraction of organic carbon)
    t = 365        # days

    # Fish
    Mfish = 1000   # g (Mass of fish)
    Cfood_common = 100 # ng/g (Concentration in food)
    IRfood = 20    # g/day (Ingestion rate)
    AFgills = 0.8  # Absorption fraction gills
    AFfood = 0.9   # Absorption fraction food
    Qgills = 100   # L/day (Gill flow rate)
    Cfish_initial_common = 10 # ng/g (Initial concentration in fish)

    # Chemical-specific parameters
    chemicals = {
        "PFOS": {
            "Cin": 2.6,         # ng/L
            "half_life_years": 91,
            "logKow": 4.0,
            "kelim": 0.069      # days^-1
        },
        "PFOA": {
            "Cin": 211300,      # ng/L
            "half_life_years": 238,
            "logKow": 4.5,
            "kelim": 0.023      # days^-1
        }
    }

    print("--- Calculation of Chemical Concentrations after 365 Days ---\n")
    print("Part 1: Water Concentration Calculation")
    print("Equation: C(t) = (Cin * Qin / (Qout + Qin * (1 + Kd * foc))) * (1 - e^(-(k + r)t))")
    print("We assume initial water concentration C0 = 0 for a pristine environment.\n")

    final_fish_concentrations = {}

    for name, params in chemicals.items():
        # --- Water Concentration Calculation ---
        print(f"--- For {name} ---")
        print(f"Calculating {name} concentration in water, C_water_{name.lower()}(365):")
        
        logKow = params["logKow"]
        logKoc = 0.81 * logKow + 0.01
        Koc = 10**logKoc
        Kd = Koc * foc # Kd = Koc * foc is a standard relationship

        half_life_days = params["half_life_years"] * 365.25
        k = math.log(2) / half_life_days
        r = Qout / V

        # The term (Qout + Qin * (1 + Kd * foc)) is used as given.
        source_term_denominator = Qout + Qin * (1 + Kd * foc)
        source_term = (params["Cin"] * Qin) / source_term_denominator
        
        loss_rate = k + r
        c_water = source_term * (1 - math.exp(-loss_rate * t))
        
        print(f"Inputs: Cin={params['Cin']} ng/L, logKow={logKow}, half-life={params['half_life_years']} years")
        print(f"Intermediate values: logKoc={logKoc:.3f}, k={k:.2e} d^-1, r={r:.2f} d^-1")
        print(f"Result: C_water_{name.lower()}(365) = {c_water:.4f} ng/L\n")

        # --- Fish Concentration Calculation ---
        print(f"Calculating {name} concentration in fish, C_fish_{name.lower()}(365):")
        
        uptake_gills = c_water * Qgills * AFgills
        uptake_food = Cfood_common * IRfood * AFfood
        
        # K is the total uptake rate constant, (ng/g)/day
        K = (uptake_gills + uptake_food) / Mfish
        kelim = params["kelim"]
        
        # Steady-state fish concentration
        c_fish_ss = K / kelim
        
        # Final fish concentration at time t
        c_fish_final = c_fish_ss + (Cfish_initial_common - c_fish_ss) * math.exp(-kelim * t)
        final_fish_concentrations[name] = c_fish_final

        print("Based on bioaccumulation model: C_fish(t) = C_fish_ss + (C_fish(0) - C_fish_ss) * e^(-kelim * t)")
        print(f"where C_fish_ss = (Uptake_gills + Uptake_food) / (Mfish * kelim)")
        print(f"Uptake from gills = {c_water:.4f} ng/L * {Qgills} L/day * {AFgills} = {uptake_gills:.2f} ng/day")
        print(f"Uptake from food = {Cfood_common} ng/g * {IRfood} g/day * {AFfood} = {uptake_food:.2f} ng/day")
        print(f"Steady-state concentration (C_fish_ss) = ({uptake_gills:.2f} + {uptake_food:.2f}) / ({Mfish} g * {kelim} d^-1) = {c_fish_ss:.4f} ng/g")
        print(f"Final concentration C_fish_{name.lower()}(365) = {c_fish_ss:.4f} + ({Cfish_initial_common} - {c_fish_ss:.4f}) * exp(-{kelim} * {t})")
        print(f"Result: C_fish_{name.lower()}(365) = {c_fish_final:.4f} ng/g\n")

    print("--------------------------------------------------\n")
    print("Part 3: Total Accumulated Concentration in Fish")
    
    c_fish_pfos = final_fish_concentrations["PFOS"]
    c_fish_pfoa = final_fish_concentrations["PFOA"]
    total_c_fish = c_fish_pfos + c_fish_pfoa
    
    print("Total POPs concentration = C_fish_pfos(365) + C_fish_pfoa(365)")
    print(f"Total POPs concentration = {c_fish_pfos:.4f} ng/g + {c_fish_pfoa:.4f} ng/g")
    print(f"Final Answer = {total_c_fish:.4f} ng/g")
    
    return total_c_fish

if __name__ == '__main__':
    final_answer = solve_pop_accumulation()
    print(f"\n<<<Total accumulated concentration: {final_answer:.4f} ng/g>>>")
