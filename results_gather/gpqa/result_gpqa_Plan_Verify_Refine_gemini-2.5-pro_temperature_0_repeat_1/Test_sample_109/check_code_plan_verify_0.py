import math

def check_correctness():
    """
    Checks the correctness of the LLM's answer by performing all calculations
    and comparing the result.
    """
    # --- Define Observational Constraints ---
    # Keck Latitude: ~ +19.8 deg N. Can see DEC > -(90 - 19.8) = -70.2 deg
    # Paranal Latitude: ~ -24.6 deg S. Can see DEC < (90 - 24.6) = +65.4 deg
    DEC_MIN = -70.2
    DEC_MAX = 65.4

    # Magnitude limit is the stricter of the two (HIRES)
    MAG_LIMIT = 16.0

    # --- Star Data ---
    stars = {
        "Star1": {"dec": -75.0, "M_V": 15.5, "d_pc": 10},
        "Star2": {"dec": 55.0, "m_V": 16.5},
        "Star3": {"dec": 48.0, "m_V": 15.5}, # Apparent magnitude is given, so E(B-V) and distance are distractors.
        "Star4": {"dec": -48.0, "M_V": 15.5, "d_pc": 10, "E_BV": 0.4},
        "Star5": {"dec": 60.0, "M_V": 16.5, "d_pc": 5}
    }

    # --- LLM's Answer ---
    llm_answer_option = 'B'
    options_map = {
        'A': {"Star2", "Star3"},
        'B': {"Star3", "Star5"},
        'C': {"Star4", "Star5"},
        'D': {"Star1", "Star4"}
    }
    llm_answer_set = options_map.get(llm_answer_option)

    if llm_answer_set is None:
        return f"The provided answer option '{llm_answer_option}' is not a valid choice."

    # --- Analysis ---
    analysis_results = {}
    observable_stars = set()

    for name, data in stars.items():
        # 1. Calculate final apparent magnitude
        final_m_V = None
        if "m_V" in data:
            final_m_V = data["m_V"]
        elif "M_V" in data and "d_pc" in data:
            # Distance modulus: m - M = 5*log10(d) - 5
            m_V_unreddened = data["M_V"] + 5 * math.log10(data["d_pc"]) - 5
            # Extinction: A_V = 3.1 * E(B-V)
            A_V = 3.1 * data.get("E_BV", 0.0)
            final_m_V = m_V_unreddened + A_V
        
        # 2. Check constraints
        dec_ok = DEC_MIN < data["dec"] < DEC_MAX
        mag_ok = final_m_V is not None and final_m_V < MAG_LIMIT
        
        is_observable = dec_ok and mag_ok
        analysis_results[name] = {
            "dec": data["dec"],
            "dec_ok": dec_ok,
            "final_m_V": final_m_V,
            "mag_ok": mag_ok,
            "is_observable": is_observable
        }
        if is_observable:
            observable_stars.add(name)

    # --- Verdict ---
    if observable_stars == llm_answer_set:
        return "Correct"
    else:
        reason = f"Incorrect. The provided answer claims {sorted(list(llm_answer_set))} are observable, but the analysis shows {sorted(list(observable_stars))} are.\n\n"
        reason += f"Here is the detailed breakdown based on the constraints (DEC between {DEC_MIN:.1f} and {DEC_MAX:.1f}; Apparent Mag < {MAG_LIMIT}):\n\n"
        for name, res in analysis_results.items():
            reason += f"--- {name} ---\n"
            reason += f"  - Declination: {res['dec']:.1f} -> Visible: {res['dec_ok']}\n"
            reason += f"  - Final Apparent Mag: {res['final_m_V']:.2f} -> Visible: {res['mag_ok']}\n"
            reason += f"  - Observable by both: {'YES' if res['is_observable'] else 'NO'}\n"
        return reason

# Run the check and print the result
print(check_correctness())