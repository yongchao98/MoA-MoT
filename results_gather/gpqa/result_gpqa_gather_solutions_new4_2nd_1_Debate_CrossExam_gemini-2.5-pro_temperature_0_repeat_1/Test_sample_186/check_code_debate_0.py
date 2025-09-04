import math

def check_correctness_of_astro_answer():
    """
    Checks the correctness of the answer to the ESPRESSO detection question.

    The function verifies two main criteria for each star:
    1. Visibility from Paranal Observatory based on declination.
    2. Brightness based on apparent magnitude, calculated if necessary.

    It then compares the calculated number of detectable stars to the
    provided answer's claim.
    """

    # --- Step 1: Define constants and criteria from the problem ---
    PARANAL_LATITUDE = -24.6  # degrees
    # A star is visible if its declination is less than 90 - |latitude|
    VISIBILITY_DEC_LIMIT = 90.0 - abs(PARANAL_LATITUDE)
    # A star is bright enough if its apparent magnitude V is <= 17.0
    BRIGHTNESS_LIMIT_V_MAG = 17.0
    # Absolute magnitude for the hypothetical stars
    HYPO_STAR_ABS_MAG = 15.0

    # --- Step 2: Define data for each star ---
    stars = [
        {'name': 'a) Canopus', 'dec': -52.7, 'type': 'known', 'v_mag': -0.74},
        {'name': 'b) Polaris', 'dec': 89.2, 'type': 'known', 'v_mag': 1.98},
        {'name': 'c) Star at 10 pc', 'dec': 0.0, 'type': 'hypothetical', 'dist_pc': 10},
        {'name': 'd) Star at 200 pc', 'dec': 0.0, 'type': 'hypothetical', 'dist_pc': 200},
        {'name': 'e) Star at 5 pc', 'dec': 0.0, 'type': 'hypothetical', 'dist_pc': 5},
        {'name': 'f) Star at 50 pc', 'dec': 0.0, 'type': 'hypothetical', 'dist_pc': 50},
    ]

    # --- Step 3: Evaluation Logic ---
    detectable_count = 0
    evaluation_log = []

    for star in stars:
        star_name = star['name']
        
        # Criterion 1: Check visibility
        is_visible = star['dec'] < VISIBILITY_DEC_LIMIT
        if not is_visible:
            evaluation_log.append(f"INFO: {star_name} is NOT detectable. Reason: Not visible from Paranal (Declination {star['dec']}° > limit {VISIBILITY_DEC_LIMIT:.1f}°).")
            continue

        # Criterion 2: Check brightness
        # First, determine the apparent magnitude (V)
        if star['type'] == 'known':
            v_mag = star['v_mag']
        else:  # For hypothetical stars, calculate V using the distance modulus formula
            # V = Mv + 5 * log10(d) - 5
            v_mag = HYPO_STAR_ABS_MAG + 5 * math.log10(star['dist_pc']) - 5
        
        is_bright_enough = v_mag <= BRIGHTNESS_LIMIT_V_MAG
        
        if is_bright_enough:
            detectable_count += 1
            evaluation_log.append(f"INFO: {star_name} is detectable. (Visible and V={v_mag:.2f} <= {BRIGHTNESS_LIMIT_V_MAG}).")
        else:
            evaluation_log.append(f"INFO: {star_name} is NOT detectable. Reason: Too faint (V={v_mag:.2f} > {BRIGHTNESS_LIMIT_V_MAG}).")

    # --- Step 4: Final Verification ---
    # The provided answer is <<<D>>>, which corresponds to a count of 3.
    expected_count = 3
    
    if detectable_count == expected_count:
        return "Correct"
    else:
        error_message = (
            f"Incorrect. The calculated number of detectable stars is {detectable_count}, "
            f"but the provided answer implies the number is {expected_count}.\n"
            "The discrepancy arises from the following evaluation:\n"
        )
        for log_entry in evaluation_log:
            error_message += f"- {log_entry}\n"
        return error_message

# Execute the check and print the result
result = check_correctness_of_astro_answer()
print(result)