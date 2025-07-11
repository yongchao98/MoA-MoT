import math

def calculate_fish_accumulation():
    """
    Calculates the accumulation of PFOS and PFOA in fish based on the provided environmental model.
    """
    # Step 1: Define all constants from the problem description.

    # Environment parameters
    V = 10000.0  # L (Volume of freshwater body)
    Qin = 900.0  # L/d (Inflow rate)
    Qout = 1600.0 # L/d (Outflow rate)
    foc = 0.001 # 0.1% (Fraction of organic carbon)
    t = 365.0  # days (Time period)
    C0 = 0.0     # ng/L (Initial concentration in water, assumed to be zero)

    # Fish parameters
    Mfish = 1000.0  # g (Average fish weight)
    Cfood = 100.0   # ng/g (Concentration of chemicals in food)
    IRfood = 20.0   # g/day (Ingestion rate of food)
    AFgills = 0.8   # (Absorption fraction through gills)
    AFfood = 0.9    # (Absorption fraction for food)
    Qgills = 100.0  # L/day (Gill flow rate)
    Cfish_given = 10.0 # ng/g (Concentration of chemicals in fish)

    # Chemical-specific parameters
    # PFOS
    params_pfos = {
        "name": "PFOS",
        "Cin": 2.6,        # ng/L (Inflow concentration)
        "t_half_years": 91.0, # years (Half-life)
        "logKow": 4.0,     # (log Octanol-water partition coefficient)
        "kelim": 0.069     # days⁻¹ (Elimination rate constant)
    }
    # PFOA
    params_pfoa = {
        "name": "PFOA",
        "Cin": 211300.0,
        "t_half_years": 238.0,
        "logKow": 4.5,
        "kelim": 0.023
    }

    results = {}

    for params in [params_pfos, params_pfoa]:
        name = params["name"]
        Cin = params["Cin"]
        t_half_years = params["t_half_years"]
        logKow = params["logKow"]
        kelim = params["kelim"]

        print(f"--- Calculations for {name} ---")

        # Step 2: Calculate derived constants needed for the water concentration formula.
        t_half_days = t_half_years * 365.0
        k = math.log(2) / t_half_days
        logKoc = 0.81 * logKow + 0.01
        Koc = 10**logKoc
        r = Qout / V

        # Step 3: Calculate the chemical concentration in the water C(t) at t=365 days.
        # As per the plan, interpreting the ambiguous 'Kd' in the given formula as 'Koc'.
        Koc_foc_term = 1 + (Koc * foc)
        denominator_Ct = Qout + Qin * Koc_foc_term
        steady_state_term = (Cin * Qin) / denominator_Ct
        exponent_term = (k + r) * t
        transient_term = 1 - math.exp(-exponent_term)
        C_water_final = steady_state_term * transient_term

        print(f"\nEquation for water concentration C(t):")
        print(f"C(t) = (Cin × Qin / (Qout + Qin × (1 + Koc × foc))) × (1−e^−(k + r)t)")
        print(f"For {name}:")
        print(f"C({t}) = ({Cin} × {Qin} / ({Qout} + {Qin} × (1 + {Koc:.2f} × {foc}))) × (1−e^−(({k:.6f} + {r:.2f}) × {t}))")
        print(f"C({t}) = {steady_state_term:.4f} ng/L × {transient_term:.4f}")
        print(f"Resulting Water Concentration for {name}: {C_water_final:.4f} ng/L")

        # Step 4: Calculate the net accumulation rate in fish.
        uptake_gills = C_water_final * Qgills * AFgills
        uptake_food = Cfood * IRfood * AFfood
        elimination = kelim * Cfish_given * Mfish
        accumulation_rate = uptake_gills + uptake_food - elimination

        print(f"\nEquation for POP accumulation rate:")
        print("Accumulation = C(t)×Qgills×AFgills + Cfood×IRfood×AFfood − kelim×Cfish×Mfish")
        print(f"For {name}:")
        print(f"Accumulation = {C_water_final:.4f}×{Qgills}×{AFgills} + {Cfood}×{IRfood}×{AFfood} − {kelim}×{Cfish_given}×{Mfish}")
        print(f"Accumulation = {uptake_gills:.2f} + {uptake_food:.2f} - {elimination:.2f}")
        print(f"Resulting Net Accumulation Rate for {name}: {accumulation_rate:.2f} ng/day\n")
        
        results[name] = accumulation_rate

    # Final summary
    print("--- FINAL SUMMARY ---")
    print(f"The net accumulation rate of PFOS in fish after 365 days is {results['PFOS']:.2f} ng/day.")
    print(f"The net accumulation rate of PFOA in fish after 365 days is {results['PFOA']:.2f} ng/day.")
    
    # Required final answer format
    # Because the question is ambiguous asking for a single "concentration accumulated", 
    # and the two results are vastly different, providing the larger, dominant value is a reasonable interpretation.
    # The PFOA accumulation rate is the dominant factor.
    final_answer = results['PFOA']
    return final_answer

if __name__ == '__main__':
    final_result = calculate_fish_accumulation()
    # The final answer is a numeric value, but printing the summary is more informative.
    # To comply with the requested format, we output the most significant calculated value.
    # print(f"<<<{final_result:.2f}>>>")
    # Let's print both, since "what are" implies multiple
    # Let's just follow the last instruction. Print and don't ask user to copy paste
    # The script has already printed everything. No need for more prints here.
    pass
    
calculate_fish_accumulation()