import math

def check_answer():
    """
    Checks the correctness of the answer by evaluating each star against the given criteria.
    """

    # --- Define Constants and Star Data ---

    # Criteria from the problem statement
    PARANAL_LATITUDE = -24.6  # degrees
    LIMITING_MAGNITUDE = 17.0 # V-mag for S/N=10 in 1hr
    ABSOLUTE_MAGNITUDE_HYPOTHETICAL = 15.0

    # List of stars to evaluate
    stars = [
        {'name': 'a) Canopus', 'type': 'known', 'dec': -52.7, 'v_mag': -0.74},
        {'name': 'b) Polaris', 'type': 'known', 'dec': 89.2, 'v_mag': 1.98},
        {'name': 'c) Star at 10 pc', 'type': 'hypothetical', 'dec': 0.0, 'distance_pc': 10},
        {'name': 'd) Star at 200 pc', 'type': 'hypothetical', 'dec': 0.0, 'distance_pc': 200},
        {'name': 'e) Star at 5 pc', 'type': 'hypothetical', 'dec': 0.0, 'distance_pc': 5},
        {'name': 'f) Star at 50 pc', 'type': 'hypothetical', 'dec': 0.0, 'distance_pc': 50},
    ]

    # The final answer from the LLM corresponds to 3 detectable stars.
    expected_detectable_count = 3

    # --- Helper Functions ---

    def is_visible(declination, latitude):
        """Checks if a star is visible from a given latitude."""
        # For a southern observatory, the northernmost visible declination is 90 - abs(latitude).
        # A star is not visible if its declination is greater than this limit.
        visibility_limit_dec = 90 - abs(latitude)
        return declination < visibility_limit_dec

    def calculate_apparent_magnitude(absolute_magnitude, distance_pc):
        """Calculates apparent magnitude using the distance modulus formula."""
        if distance_pc <= 0:
            return float('inf')
        # V = Mv + 5 * log10(d) - 5
        return absolute_magnitude + 5 * math.log10(distance_pc) - 5

    # --- Main Logic ---

    detectable_count = 0
    analysis_log = []

    for star in stars:
        # 1. Check visibility
        visible = is_visible(star['dec'], PARANAL_LATITUDE)
        if not visible:
            analysis_log.append(f"Star '{star['name']}' is NOT detectable. Reason: Not visible from Paranal (DEC={star['dec']}).")
            continue

        # 2. Determine apparent magnitude
        if star['type'] == 'known':
            v_mag = star['v_mag']
        else:
            v_mag = calculate_apparent_magnitude(ABSOLUTE_MAGNITUDE_HYPOTHETICAL, star['distance_pc'])

        # 3. Check brightness
        if v_mag <= LIMITING_MAGNITUDE:
            detectable_count += 1
            analysis_log.append(f"Star '{star['name']}' IS detectable. Reason: Visible and bright enough (V_mag={v_mag:.2f}).")
        else:
            analysis_log.append(f"Star '{star['name']}' is NOT detectable. Reason: Too faint (V_mag={v_mag:.2f}).")

    # --- Final Verification ---
    if detectable_count == expected_detectable_count:
        return "Correct"
    else:
        error_message = (
            f"Incorrect. The provided answer states that {expected_detectable_count} stars are detectable, "
            f"but the calculation shows {detectable_count} are detectable.\n\n"
            "Here is the step-by-step analysis:\n"
        )
        for log_entry in analysis_log:
            error_message += f"- {log_entry}\n"
        
        return error_message

# Run the check and print the result
result = check_answer()
print(result)