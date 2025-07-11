import math

def calculate_accumulation():
    """
    Calculates the net accumulation rate of PFOS and PFOA in fish
    based on the provided environmental and physiological parameters.
    """
    # --- General Environmental and Fish Parameters ---
    V_water = 10000  # L
    Q_in = 900  # L/d
    Q_out = 1600  # L/d
    foc = 0.001  # 0.1% organic carbon
    t = 365  # days
    C0 = 0 # ng/L, assuming pristine environment initially

    M_fish = 1000  # g
    IR_food = 20  # g/day
    AF_gills = 0.8
    AF_food = 0.9
    Q_gills = 100  # L/day

    # --- Chemical-Specific Parameters ---
    # PFOS
    C_in_pfos = 2.6  # ng/L
    half_life_pfos_years = 91  # years
    log_Kow_pfos = 4.0
    C_food_pfos = 100  # ng/g
    k_elim_pfos = 0.069  # days^-1
    C_fish_pfos = 10  # ng/g

    # PFOA
    C_in_pfoa = 211300  # ng/L
    half_life_pfoa_years = 238  # years
    log_Kow_pfoa = 4.5
    C_food_pfoa = 100  # ng/g
    k_elim_pfoa = 0.023  # days^-1
    C_fish_pfoa = 10  # ng/g

    chemicals = [
        {
            "name": "PFOS",
            "C_in": C_in_pfos,
            "half_life_years": half_life_pfos_years,
            "log_Kow": log_Kow_pfos,
            "C_food": C_food_pfos,
            "k_elim": k_elim_pfos,
            "C_fish": C_fish_pfos
        },
        {
            "name": "PFOA",
            "C_in": C_in_pfoa,
            "half_life_years": half_life_pfoa_years,
            "log_Kow": log_Kow_pfoa,
            "C_food": C_food_pfoa,
            "k_elim": k_elim_pfoa,
            "C_fish": C_fish_pfoa
        }
    ]

    results = {}

    for chem in chemicals:
        print(f"--- Calculating for {chem['name']} ---")

        # Step 1: Calculate intermediate values for C(t)
        log_Koc = 0.81 * chem['log_Kow'] + 0.01
        Koc = 10**log_Koc
        Kd = Koc * foc
        half_life_days = chem['half_life_years'] * 365.25
        k = math.log(2) / half_life_days
        r = Q_out / V_water
        
        # Step 2: Calculate C(t)
        # C(t) = (Cin * Qin / (Qout + Qin * (1 + Kd * foc))) * (1 - e^-(k + r)t) + C0 * e^-(k + r)t
        exponent = -(k + r) * t
        exp_term = math.exp(exponent)
        
        # Since C0 is 0, the second part of the equation is 0.
        # The exp_term is very close to 0 for t=365, so (1 - exp_term) is ~1.
        c_t_numerator = chem['C_in'] * Q_in
        c_t_denominator = (Q_out + Q_in * (1 + Kd * foc))
        C_t = (c_t_numerator / c_t_denominator) * (1 - exp_term)

        print(f"Water concentration C(t) for {chem['name']} at {t} days:")
        print(f"C({t}) = ({chem['C_in']} * {Q_in} / ({Q_out} + {Q_in} * (1 + {Kd:.4f} * {foc}))) * (1 - e^-(({k:.6f} + {r}) * {t}))")
        print(f"C({t}) = {C_t:.4f} ng/L\n")

        # Step 3: Calculate POP accumulation rate
        # POP accumulation = C(t) * Qgills * AFgills + Cfood * IRfood * AFfood - kelim * Cfish * Mfish
        gill_uptake = C_t * Q_gills * AF_gills
        food_uptake = chem['C_food'] * IR_food * AF_food
        elimination = chem['k_elim'] * chem['C_fish'] * M_fish
        
        net_accumulation_rate = gill_uptake + food_uptake - elimination
        results[chem['name']] = net_accumulation_rate

        print(f"Net accumulation rate for {chem['name']} in fish at {t} days:")
        print(f"Accumulation Rate = C(t) * Qgills * AFgills + Cfood * IRfood * AFfood - kelim * Cfish * Mfish")
        print(f"Accumulation Rate = {C_t:.4f} * {Q_gills} * {AF_gills} + {chem['C_food']} * {IR_food} * {AF_food} - {chem['k_elim']} * {chem['C_fish']} * {M_fish}")
        print(f"Accumulation Rate = {gill_uptake:.2f} + {food_uptake:.2f} - {elimination:.2f}")
        print(f"Net Accumulation Rate for {chem['name']}: {net_accumulation_rate:.2f} ng/day\n")

    return results

# Run the calculation and print the final answer in the specified format
final_results = calculate_accumulation()
pfos_result = final_results.get("PFOS", 0)
pfoa_result = final_results.get("PFOA", 0)

print(f"<<<PFOS accumulation rate: {pfos_result:.2f} ng/day, PFOA accumulation rate: {pfoa_result:.2f} ng/day>>>")
