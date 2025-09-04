import math

def check_star_observability():
    """
    Checks which stars are observable by both ESPRESSO and HIRES based on
    declination and apparent magnitude constraints.
    """
    # --- 1. Define Constraints ---
    # Magnitude limit is the stricter of the two (V < 17 for ESPRESSO, V < 16 for HIRES)
    MAG_LIMIT = 16.0

    # Declination limits based on observatory latitudes
    # Keck Observatory (HIRES) Latitude: ~+19.8 deg N
    # Paranal Observatory (ESPRESSO) Latitude: ~-24.6 deg S
    # A star is visible if its declination is between (lat - 90) and (lat + 90)
    DEC_MIN = 19.8 - 90.0  # Keck's southern limit
    DEC_MAX = -24.6 + 90.0 # Paranal's northern limit

    # --- 2. Define Star Data ---
    # Note: RA in hours is converted to degrees for consistency, though not used in the check.
    stars = {
        "Star1": {"DEC": -75.0, "M_V": 15.5, "dist_pc": 10, "E_BV": 0.0},
        "Star2": {"DEC": 55.0, "V_app": 16.5},
        "Star3": {"DEC": 48.0, "V_app": 15.5}, # Assuming V_app is the final observed magnitude
        "Star4": {"DEC": -48.0, "M_V": 15.5, "dist_pc": 10, "E_BV": 0.4},
        "Star5": {"DEC": 60.0, "M_V": 16.5, "dist_pc": 5, "E_BV": 0.0},
    }

    # --- 3. Analyze Each Star ---
    observable_stars = []
    analysis_log = []

    for name, data in stars.items():
        # Check Declination
        dec = data["DEC"]
        dec_visible = DEC_MIN < dec < DEC_MAX

        # Check Apparent Magnitude
        if "V_app" in data:
            v_mag = data["V_app"]
        else:
            # Calculate apparent magnitude: V = M_V + 5*log10(d) - 5 + A_V
            # A_V = 3.1 * E(B-V)
            a_v = 3.1 * data["E_BV"]
            v_mag = data["M_V"] + 5 * math.log10(data["dist_pc"]) - 5 + a_v
        
        mag_observable = v_mag < MAG_LIMIT

        # Store results
        if dec_visible and mag_observable:
            observable_stars.append(name)
        
        # Log the reasoning for each star
        log_entry = f"{name}: DEC={dec:.1f}, V_app={v_mag:.2f}. "
        if not dec_visible:
            log_entry += f"FAIL (Declination not in range {DEC_MIN:.1f} to {DEC_MAX:.1f})."
        elif not mag_observable:
            log_entry += f"FAIL (Magnitude not < {MAG_LIMIT})."
        else:
            log_entry += "PASS."
        analysis_log.append(log_entry)

    # --- 4. Verify the LLM's Answer ---
    llm_answer_choice = "D"
    llm_answer_stars = ["Star3", "Star5"]

    # Sort for consistent comparison
    observable_stars.sort()
    llm_answer_stars.sort()

    if observable_stars == llm_answer_stars:
        return "Correct"
    else:
        error_message = f"Incorrect. The LLM's answer is {llm_answer_choice} ({', '.join(llm_answer_stars)}), but the correct set of observable stars is ({', '.join(observable_stars)}).\n\n"
        error_message += "Detailed Analysis:\n"
        error_message += "\n".join(analysis_log)
        return error_message

# Run the check and print the result
result = check_star_observability()
print(result)