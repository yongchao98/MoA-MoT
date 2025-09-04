import math

def check_espresso_detectability():
    """
    Checks how many stars from a given list are detectable by ESPRESSO at Paranal,
    and verifies if the LLM's answer is correct.
    """
    # --- Define Constraints and Constants ---
    # Constraint 1: Observability from Paranal
    PARANAL_LATITUDE = -24.6  # degrees
    # A star is observable if its declination is less than 90 - |latitude|
    MAX_OBSERVABLE_DEC = 90 - abs(PARANAL_LATITUDE)

    # Constraint 2: Brightness for S/N >= 10 in 1 hour
    # From ESPRESSO performance table, V=15 gives S/N=11, V=16 gives S/N=7.
    # Interpolating for S/N=10: m_lim = 15 + (10-11)/(7-11) = 15.25
    LIMITING_MAGNITUDE = 15.25

    # The LLM's answer is C, which corresponds to 3 stars.
    llm_answer_value = 3

    # --- Star Data ---
    stars = {
        'a) Canopus': {'dec': -52.7, 'm_V': -0.74, 'M_V': None, 'dist_pc': None},
        'b) Polaris': {'dec': 89.26, 'm_V': 1.98, 'M_V': None, 'dist_pc': None},
        'c) Star at 10 pc': {'dec': 0.0, 'm_V': None, 'M_V': 15.0, 'dist_pc': 10},
        'd) Star at 200 pc': {'dec': 0.0, 'm_V': None, 'M_V': 15.0, 'dist_pc': 200},
        'e) Star at 5 pc': {'dec': 0.0, 'm_V': None, 'M_V': 15.0, 'dist_pc': 5},
        'f) Star at 50 pc': {'dec': 0.0, 'm_V': None, 'M_V': 15.0, 'dist_pc': 50},
    }

    detectable_count = 0
    evaluation_log = []
    detectable_stars = []

    for name, data in stars.items():
        # --- Constraint 1: Observability ---
        is_observable = data['dec'] <= MAX_OBSERVABLE_DEC
        if not is_observable:
            evaluation_log.append(f"FAIL: {name} is not observable (DEC {data['dec']:.2f}° > {MAX_OBSERVABLE_DEC:.2f}°).")
            continue

        # --- Constraint 2: Brightness ---
        apparent_magnitude = data['m_V']
        # Calculate apparent magnitude if not given
        if apparent_magnitude is None:
            # m_V = M_V + 5 * log10(d) - 5
            apparent_magnitude = data['M_V'] + 5 * math.log10(data['dist_pc']) - 5

        is_bright_enough = apparent_magnitude <= LIMITING_MAGNITUDE
        if not is_bright_enough:
            evaluation_log.append(f"FAIL: {name} is too faint (m_V {apparent_magnitude:.2f} > {LIMITING_MAGNITUDE:.2f}).")
            continue
        
        # If both constraints pass, the star is detectable
        evaluation_log.append(f"PASS: {name} is detectable (Observable and m_V = {apparent_magnitude:.2f}).")
        detectable_count += 1
        detectable_stars.append(name)

    # --- Final Verification ---
    if detectable_count == llm_answer_value:
        return "Correct"
    else:
        reason = f"The answer is incorrect. The calculated number of detectable stars is {detectable_count}, but the provided answer corresponds to {llm_answer_value}.\n"
        reason += f"The stars calculated to be detectable are: {detectable_stars}\n"
        reason += "Detailed evaluation log:\n"
        for log_entry in evaluation_log:
            reason += f"- {log_entry}\n"
        return reason

# Execute the check and print the result
result = check_espresso_detectability()
print(result)