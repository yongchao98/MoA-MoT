import math

def check_correctness():
    """
    Checks the correctness of the answer by verifying the number of detectable stars.
    A star is detectable if it is:
    1. Visible from Paranal Observatory (DEC < +65.4 deg).
    2. Bright enough for ESPRESSO (Apparent Magnitude V <= 17.0).
    """

    # --- Define Constraints ---
    # Visibility constraint for Paranal Observatory (latitude ~24.6 S)
    paranal_latitude = -24.6
    visibility_limit_dec = 90 - abs(paranal_latitude)

    # Brightness constraint based on ESPRESSO performance for S/N >= 10
    # V=16 -> S/N=15; V=19 -> S/N=3. The limit is ~V=17.0.
    limiting_magnitude = 17.0

    # --- Define Star Data ---
    stars = [
        {'id': 'a) Canopus', 'dec': -52.7, 'V': -0.74, 'M_V': None, 'd_pc': None},
        {'id': 'b) Polaris', 'dec': 89.3, 'V': 1.98, 'M_V': None, 'd_pc': None},
        {'id': 'c) Star at 10 pc', 'dec': 0, 'V': None, 'M_V': 15, 'd_pc': 10},
        {'id': 'd) Star at 200 pc', 'dec': 0, 'V': None, 'M_V': 15, 'd_pc': 200},
        {'id': 'e) Star at 5 pc', 'dec': 0, 'V': None, 'M_V': 15, 'd_pc': 5},
        {'id': 'f) Star at 50 pc', 'dec': 0, 'V': None, 'M_V': 15, 'd_pc': 50},
    ]

    detectable_stars_ids = []
    analysis_log = []

    # --- Main Analysis Loop ---
    for star in stars:
        # 1. Check Visibility
        is_visible = star['dec'] <= visibility_limit_dec
        if not is_visible:
            analysis_log.append(f"{star['id']}: Not detectable. Reason: Not visible from Paranal (DEC={star['dec']} > {visibility_limit_dec:.1f}).")
            continue

        # 2. Check Brightness
        # Calculate apparent magnitude if not provided
        apparent_magnitude = star['V']
        if apparent_magnitude is None:
            M_V = star['M_V']
            d_pc = star['d_pc']
            # Distance modulus formula: V = M + 5*log10(d/10)
            apparent_magnitude = M_V + 5 * math.log10(d_pc / 10)

        is_bright_enough = apparent_magnitude <= limiting_magnitude
        if not is_bright_enough:
            analysis_log.append(f"{star['id']}: Not detectable. Reason: Too faint (V={apparent_magnitude:.2f} > {limiting_magnitude}).")
            continue
        
        # If both checks pass, the star is detectable
        analysis_log.append(f"{star['id']}: Detectable. (Visible and V={apparent_magnitude:.2f} <= {limiting_magnitude}).")
        detectable_stars_ids.append(star['id'])

    # --- Final Verification ---
    # The provided answer is 'D', which corresponds to 3.
    expected_count = 3
    
    if len(detectable_stars_ids) == expected_count:
        # Further check if the correct stars were identified
        expected_stars = ['a) Canopus', 'c) Star at 10 pc', 'e) Star at 5 pc']
        if sorted(detectable_stars_ids) == sorted(expected_stars):
            return "Correct"
        else:
            return (f"Incorrect. The number of detectable stars is correct ({expected_count}), "
                    f"but the specific stars identified are wrong. "
                    f"Expected: {sorted(expected_stars)}, Found: {sorted(detectable_stars_ids)}.")
    else:
        return (f"Incorrect. The final answer claims {expected_count} stars are detectable, "
                f"but the calculation shows {len(detectable_stars_ids)} are. "
                f"The detectable stars found were: {detectable_stars_ids}.")

# Run the check
result = check_correctness()
print(result)