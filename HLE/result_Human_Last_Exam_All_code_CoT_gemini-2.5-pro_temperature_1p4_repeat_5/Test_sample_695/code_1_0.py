import math

def calculate_fish_accumulation():
    """
    Calculates and prints the net accumulation rate of PFOS and PFOA in fish
    based on the environmental and biological parameters provided.
    """
    # --- General and Fish Parameters ---
    V_water = 10000  # L (Volume of water body)
    Qin = 900  # L/d (Inflow rate)
    Qout = 1600  # L/d (Outflow rate)
    foc = 0.001  # 0.1% (Fraction of organic carbon)
    t = 365  # days
    C0_water = 0 # ng/L (Assuming pristine environment initially)

    Mfish = 1000  # g (Fish weight)
    Cfood = 100  # ng/g (Concentration in food)
    IRfood = 20  # g/day (Ingestion rate)
    AFfood = 0.9  # (Absorption fraction, food)
    Qgills = 100  # L/day (Gill flow rate)
    AFgills = 0.8  # (Absorption fraction, gills)
    Cfish = 10  # ng/g (Initial concentration in fish)

    # --- Chemical-specific Parameters ---
    chemicals = {
        "PFOS": {
            "Cin_water": 2.6,  # ng/L
            "half_life_years": 91,
            "logKow": 4.0,
            "kelim": 0.069  # days⁻¹
        },
        "PFOA": {
            "Cin_water": 211300,  # ng/L
            "half_life_years": 238,
            "logKow": 4.5,
            "kelim": 0.023  # days⁻¹
        }
    }

    # --- Calculations for each chemical ---
    results = {}
    for name, params in chemicals.items():
        print(f"--- Calculating accumulation for {name} ---")

        # 1. Calculate water concentration C(t)
        # log Koc = 0.81 * log Kow + 0.01
        logKoc = 0.81 * params["logKow"] + 0.01
        Koc = 10**logKoc
        # Kd = Koc * foc
        Kd = Koc * foc
        
        # k = ln(2) / half_life_in_days
        half_life_days = params["half_life_years"] * 365
        k = math.log(2) / half_life_days
        
        # r = Qout / V
        r = Qout / V_water
        
        k_plus_r = k + r

        # C(t) = (Cin * Qin / (Qout + Qin * (1 + Kd * foc))) * (1 - e^-(k + r)t) + C0 * e^-(k + r)t
        # Since C0 = 0, the second term is zero.
        water_conc_steady_state_part_numerator = params["Cin_water"] * Qin
        water_conc_steady_state_part_denominator = Qout + Qin * (1 + Kd * foc)
        A = water_conc_steady_state_part_numerator / water_conc_steady_state_part_denominator
        
        exp_term = math.exp(-k_plus_r * t)
        C_t_water = A * (1 - exp_term)

        print("1. Water Concentration C(t) after 365 days:")
        print(f"   logKoc = 0.81 * {params['logKow']} + 0.01 = {logKoc:.3f}")
        print(f"   Koc = 10^{logKoc:.3f} = {Koc:.2f}")
        print(f"   Kd = {Koc:.2f} * {foc} = {Kd:.3f}")
        print(f"   k = ln(2) / ({params['half_life_years']} * 365) = {k:.3e} day⁻¹")
        print(f"   r = {Qout} L/d / {V_water} L = {r:.3f} day⁻¹")
        print("   Equation: C(t) = (Cin * Qin / (Qout + Qin * (1 + Kd * foc))) * (1 - e^-(k + r)t)")
        print(f"   C(365) = ({params['Cin_water']} * {Qin} / ({Qout} + {Qin} * (1 + {Kd:.3f} * {foc}))) * (1 - e^-({k:.3e} + {r:.3f})*{t})")
        print(f"   C(365) = {C_t_water:.3f} ng/L\n")

        # 2. Calculate POP accumulation in fish
        uptake_gills = C_t_water * Qgills * AFgills
        uptake_food = Cfood * IRfood * AFfood
        elimination = params["kelim"] * Cfish * Mfish
        net_accumulation = uptake_gills + uptake_food - elimination
        results[name] = net_accumulation

        print("2. Net Accumulation Rate in Fish:")
        print("   Equation: Accumulation = C(water) * Qgills * AFgills + Cfood * IRfood * AFfood - kelim * Cfish * Mfish")
        print(f"   Uptake (Gills) = {C_t_water:.3f} ng/L * {Qgills} L/day * {AFgills} = {uptake_gills:.2f} ng/day")
        print(f"   Uptake (Food)  = {Cfood} ng/g * {IRfood} g/day * {AFfood} = {uptake_food:.2f} ng/day")
        print(f"   Elimination    = {params['kelim']} day⁻¹ * {Cfish} ng/g * {Mfish} g = {elimination:.2f} ng/day")
        print(f"   Net Accumulation ({name}) = {uptake_gills:.2f} + {uptake_food:.2f} - {elimination:.2f} = {net_accumulation:.2f} ng/day\n")

    pfos_result = results.get("PFOS", 0)
    pfoa_result = results.get("PFOA", 0)
    
    final_answer = f"The net accumulation rate for PFOS is {pfos_result:.2f} ng/day, and for PFOA is {pfoa_result:.2f} ng/day."
    # The problem format requires a specific end tag.
    # While the calculation provides two distinct results, the final format expects a single content block.
    # We will combine the two results into one summary string.
    print("--- Final Answer ---")
    print(final_answer)
    return final_answer

# Execute the calculation and print the final answer in the required format
final_answer_string = calculate_fish_accumulation()
print(f"<<<{final_answer_string}>>>")