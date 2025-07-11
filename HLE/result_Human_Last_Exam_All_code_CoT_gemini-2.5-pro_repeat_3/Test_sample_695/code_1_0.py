import math

def calculate_fish_accumulation():
    """
    Calculates the accumulation rate of PFOS and PFOA in fish after 365 days.
    """
    # --- 1. Define Constants ---
    # Environment
    V = 10000  # L (Volume of water body)
    Qin = 900  # L/d (Inflow rate)
    Qout = 1600 # L/d (Outflow rate)
    foc = 0.001 # 0.1% organic carbon
    C0 = 0      # ng/L (Initial concentration in water, assumed to be 0)
    t = 365     # days

    # Fish
    Mfish = 1000   # g (Fish weight)
    Cfood = 100    # ng/g (Food concentration for each chemical)
    IRfood = 20    # g/day (Ingestion rate)
    AFgills = 0.8  # Absorption fraction, gills
    AFfood = 0.9   # Absorption fraction, food
    Qgills = 100   # L/day (Gill flow rate)
    Cfish_initial = 10 # ng/g (Initial fish concentration for both chemicals)

    # Chemical-specific parameters
    params = {
        "PFOS": {
            "Cin": 2.6,        # ng/L
            "half_life_y": 91, # years
            "logKow": 4.0,
            "kelim": 0.069     # days⁻¹
        },
        "PFOA": {
            "Cin": 211300,     # ng/L
            "half_life_y": 238,# years
            "logKow": 4.5,
            "kelim": 0.023     # days⁻¹
        }
    }

    print("This script calculates the chemical accumulation rate in fish after 365 days.")
    print("---")

    results = {}

    for chemical, p in params.items():
        print(f"Calculating for {chemical}:")

        # --- 2. Calculate Water Concentration C(t) ---
        # Convert half-life from years to days
        half_life_d = p["half_life_y"] * 365.0
        # Calculate degradation rate constant k
        k = math.log(2) / half_life_d
        # Calculate log Koc
        logKoc = 0.81 * p["logKow"] + 0.01
        # Calculate Koc
        Koc = 10**logKoc
        # Calculate Kd
        Kd = Koc * foc
        # Calculate hydraulic removal rate r
        r = Qout / V
        # Total removal rate for exponent
        k_total = k + r

        # Calculate the steady-state factor from the provided equation
        # Denominator = Qout + Qin * (1 + Kd * foc)
        # Note: The units in this part of the user-provided formula are inconsistent.
        # We proceed using the formula as given.
        denominator = Qout + Qin * (1 + Kd * foc)
        c_ss_factor = (p["Cin"] * Qin) / denominator

        # Calculate exponent term e^-(k+r)t
        exp_term = math.exp(-k_total * t)

        # Calculate water concentration at t=365 days
        C_t = c_ss_factor * (1 - exp_term) + C0 * exp_term

        # --- 3. Calculate Chemical Accumulation Rate in Fish ---
        # Uptake from gills (ng/day)
        uptake_gills = C_t * Qgills * AFgills
        # Uptake from food (ng/day)
        uptake_food = Cfood * IRfood * AFfood
        # Elimination (ng/day)
        elimination = p["kelim"] * Cfish_initial * Mfish
        # Net accumulation rate (ng/day)
        accumulation_rate = uptake_gills + uptake_food - elimination
        results[chemical] = accumulation_rate

        # --- 4. Display Final Equations and Results ---
        print(f"Step 1: Water Concentration of {chemical} at {t} days, C({t})")
        print(f"C({t}) = ({p['Cin']} * {Qin} / ({Qout} + {Qin} * (1 + {Koc:.3f} * {foc}))) * (1 - e^-(({k:.5e} + {r}) * {t})) + {C0} * e^-(({k:.5e} + {r}) * {t})")
        print(f"C({t}) = {C_t:.3f} ng/L\n")

        print(f"Step 2: Accumulation Rate of {chemical} in fish at {t} days")
        print(f"Rate = C(t) * Qgills * AFgills + Cfood * IRfood * AFfood - kelim * Cfish * Mfish")
        print(f"Rate = {C_t:.3f} * {Qgills} * {AFgills} + {Cfood} * {IRfood} * {AFfood} - {p['kelim']} * {Cfish_initial} * {Mfish}")
        print(f"Rate = {uptake_gills:.1f} + {uptake_food:.1f} - {elimination:.1f}")
        print(f"The net accumulation rate for {chemical} is {accumulation_rate:.1f} ng/day.\n---")

    print(f"Final Answer:")
    print(f"The accumulation rate for PFOS at 365 days is {results['PFOS']:.1f} ng/day.")
    print(f"The accumulation rate for PFOA at 365 days is {results['PFOA']:.1f} ng/day.")


if __name__ == '__main__':
    calculate_fish_accumulation()
    # The final answer is the set of calculated accumulation rates.
    # PFOS: 1184.8 ng/day
    # PFOA: 6077130.4 ng/day
    # The question asks for a single value "concentration accumulated", which is ambiguous.
    # Returning the sum of the rates as a potential interpretation.
    # Or maybe the larger of the two. Let's return the two distinct values as the most complete answer.
    # Since the prompt asks for a single return value in the format <<<answer>>>,
    # and the question uses "chemical concentration" (singular), I'll return the value for PFOA
    # as it is the dominant contributor by several orders of magnitude.
    # But the prompt also says "what are the chemical concentration" (plural).
    # This is contradictory.
    # The best approach is to provide both numbers in the final output text,
    # and then perhaps a combined or most significant value in the <<<>>> format.
    # Given the ambiguity, I'll provide both results clearly in the text output.
    # The problem asks for "the chemical concentration accumulated". I will provide both values as the answer.
    # Let's provide the PFOA value as it's the more significant one.
    # No, let's provide both in a string.
    
    final_answer = "PFOS Rate: 1184.8 ng/day, PFOA Rate: 6077130.4 ng/day"
    # The example format is <<<C>>> or <<<9.8>>>, suggesting a single letter or number.
    # Given the ambiguity, I will provide the larger number.
    # PFOA accumulation rate is 6077130.4 ng/day.
    final_pfoa_rate = 6077130.4
    # The user might be expecting the total accumulated, not a rate.
    # But the formula is a rate. Sticking with the rate.
    # print(f"<<<{final_pfoa_rate:.1f}>>>")
    # Let's just output the text as requested. The <<<>>> is a final step.
    # PFOS accumulation rate is 1184.83 ng/day
    # PFOA accumulation rate is 6077130.4 ng/day
    # I will present both answers clearly.