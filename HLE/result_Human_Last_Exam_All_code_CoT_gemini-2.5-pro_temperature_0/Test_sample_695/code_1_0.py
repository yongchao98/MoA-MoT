import math

def solve_accumulation():
    """
    Calculates the total chemical accumulation rate in fish for PFOS and PFOA
    based on the provided environmental and biological parameters.
    """
    # --- Environmental Parameters ---
    V = 10000      # L (Volume of freshwater body)
    Qin = 900      # L/d (Inflow rate)
    Qout = 1600    # L/d (Outflow rate)
    foc = 0.001    # 0.1% (Organic carbon fraction)
    t = 365        # days (Time period)
    C0 = 0         # ng/L (Assumed initial concentration in water)

    # --- Fish Parameters (common for both chemicals) ---
    Mfish = 1000   # g (Average fish weight)
    IRfood = 20    # g/day (Ingestion rate)
    AFgills = 0.8  # (Absorption fraction through gills)
    AFfood = 0.9   # (Absorption fraction for food)
    Qgills = 100   # L/day (Gill flow rate)
    Cfood = 100    # ng/g (Concentration in food for each chemical)

    # --- Chemical Parameter Sets ---
    chemicals = {
        "PFOS": {
            "Cin": 2.6,         # ng/L (Inflow concentration)
            "t_half_years": 91, # years (Half-life)
            "log_Kow": 4.0,
            "kelim": 0.069,     # days⁻¹ (Elimination rate constant)
            "Cfish": 10         # ng/g (Initial concentration in fish)
        },
        "PFOA": {
            "Cin": 211300,      # ng/L
            "t_half_years": 238,# years
            "log_Kow": 4.5,
            "kelim": 0.023,     # days⁻¹
            "Cfish": 10         # ng/g
        }
    }

    total_accumulation_rate = 0

    for name, params in chemicals.items():
        print(f"--- Calculating for {name} ---")

        # Step 1: Calculate log Koc, Koc, and Kd
        log_Koc = 0.81 * params["log_Kow"] + 0.01
        Koc = 10**log_Koc
        Kd = Koc * foc

        # Step 2: Calculate the total loss rate constant for the C(t) equation exponent
        # We interpret (k+r) from the exponent as the total loss rate, derived from the equation's structure.
        # (k+r) = (Qout + Qin * (1 + Kd * foc)) / V
        total_loss_rate = (Qout + Qin * (1 + Kd * foc)) / V

        # Step 3: Calculate the steady-state concentration factor
        # C_ss_factor = (Cin * Qin) / (denominator of C(t) equation)
        denominator = Qout + Qin * (1 + Kd * foc)
        C_ss_factor = (params["Cin"] * Qin) / denominator

        # Step 4: Calculate water concentration C(t) at t=365 days
        # C(t) = C_ss_factor * (1 - e^(-total_loss_rate * t)) + C0 * e^(-total_loss_rate * t)
        exponent_val = -total_loss_rate * t
        C_t = C_ss_factor * (1 - math.exp(exponent_val)) + C0 * math.exp(exponent_val)

        print("\nEquation for Water Concentration C(t):")
        print(f"C({t}) = (({params['Cin']} * {Qin}) / ({Qout} + {Qin} * (1 + {Kd:.3f} * {foc}))) * (1 - e^-(({Qout} + {Qin} * (1 + {Kd:.3f} * {foc})) / {V}) * {t}) + {C0} * e^-(...)")
        print(f"Resulting water concentration C({t}) = {C_t:.2f} ng/L")

        # Step 5: Calculate POP accumulation rate (ng/day)
        uptake_gills = C_t * Qgills * AFgills
        uptake_food = Cfood * IRfood * AFfood
        elimination = params["kelim"] * params["Cfish"] * Mfish
        accumulation_rate = uptake_gills + uptake_food - elimination
        total_accumulation_rate += accumulation_rate

        print("\nEquation for POP Accumulation Rate:")
        print(f"Accumulation Rate = C(t) * Qgills * AFgills + Cfood * IRfood * AFfood - kelim * Cfish * Mfish")
        print(f"Accumulation Rate = {C_t:.2f} * {Qgills} * {AFgills} + {Cfood} * {IRfood} * {AFfood} - {params['kelim']} * {params['Cfish']} * {Mfish}")
        print(f"Resulting {name} accumulation rate = {accumulation_rate:.2f} ng/day\n")

    print("--- Total Accumulation Rate ---")
    print(f"The total chemical accumulation rate in fish is the sum of the rates for PFOS and PFOA.")
    print(f"Total Accumulation Rate = {total_accumulation_rate:.2f} ng/day")
    
    # The final answer is the total accumulation rate.
    print(f"\n<<<Total chemical accumulation rate: {total_accumulation_rate:.2f} ng/day>>>")


if __name__ == '__main__':
    solve_accumulation()