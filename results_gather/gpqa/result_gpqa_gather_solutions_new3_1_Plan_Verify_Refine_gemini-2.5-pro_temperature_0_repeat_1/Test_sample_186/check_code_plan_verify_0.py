import math

def check_correctness_of_answer():
    """
    This function checks the correctness of the provided answer by verifying the two main constraints:
    1. Visibility of the star from Paranal Observatory.
    2. Brightness of the star compared to ESPRESSO's detection limit.
    """

    # --- Step 1: Define constants and criteria from the provided answer ---
    # Latitude of Paranal Observatory in degrees.
    PARANAL_LATITUDE = -24.6
    # The limiting apparent magnitude for S/N=10 in a 1-hour exposure, as per the solution.
    LIMITING_MAGNITUDE = 17.1
    # The expected number of detectable stars from the answer. The answer is 'C', which corresponds to 3.
    EXPECTED_COUNT = 3

    # --- Step 2: Define the star data ---
    stars = [
        {'name': 'a) Canopus', 'V': -0.74, 'DEC': -52.7},
        {'name': 'b) Polaris', 'V': 1.98, 'DEC': 89.26},
        {'name': 'c) Star at 10 pc', 'Mv': 15, 'd': 10, 'DEC': 0},
        {'name': 'd) Star at 200 pc', 'Mv': 15, 'd': 200, 'DEC': 0},
        {'name': 'e) Star at 5 pc', 'Mv': 15, 'd': 5, 'DEC': 0},
        {'name': 'f) Star at 50 pc', 'Mv': 15, 'd': 50, 'DEC': 0},
    ]

    # --- Step 3: Define helper functions ---
    def is_visible(declination, latitude):
        """Checks if a star is ever visible from a given latitude."""
        # For a southern hemisphere observatory, the northern declination limit is 90 - |latitude|.
        # A star is not visible if its declination is greater than this limit.
        visibility_limit_dec = 90 - abs(latitude)
        return declination <= visibility_limit_dec

    def calculate_apparent_magnitude(Mv, d):
        """Calculates apparent magnitude (m) from absolute magnitude (M) and distance (d) in parsecs."""
        # Formula: m = M - 5 + 5 * log10(d)
        if d <= 0:
            return float('inf')  # Avoid math domain error for non-physical distances
        return Mv - 5 + 5 * math.log10(d)

    # --- Step 4: Iterate and check each star ---
    detectable_count = 0
    detectable_stars_list = []
    reasons_for_rejection = []

    for star in stars:
        # Check 1: Visibility
        visible = is_visible(star['DEC'], PARANAL_LATITUDE)
        if not visible:
            reasons_for_rejection.append(
                f"Star '{star['name']}' is not visible from Paranal (DEC={star['DEC']:.2f}°)."
            )
            continue

        # Check 2: Brightness
        if 'V' in star:
            apparent_magnitude = star['V']
        else:
            apparent_magnitude = calculate_apparent_magnitude(star['Mv'], star['d'])

        bright_enough = apparent_magnitude <= LIMITING_MAGNITUDE
        if not bright_enough:
            reasons_for_rejection.append(
                f"Star '{star['name']}' is too faint (V={apparent_magnitude:.2f} > limit of {LIMITING_MAGNITUDE})."
            )
            continue

        # If both checks pass, the star is detectable
        detectable_count += 1
        detectable_stars_list.append(star['name'])

    # --- Step 5: Final verification ---
    if detectable_count != EXPECTED_COUNT:
        return (f"Incorrect. The answer claims {EXPECTED_COUNT} stars are detectable, but the code found {detectable_count}. "
                f"The detectable stars are: {detectable_stars_list}. "
                f"Reasons for other stars being rejected: {reasons_for_rejection}")

    # Additionally, check if the specific stars match the answer's reasoning
    expected_detectable_stars = {'a) Canopus', 'c) Star at 10 pc', 'e) Star at 5 pc'}
    if set(detectable_stars_list) != expected_detectable_stars:
        return (f"Incorrect. The number of detectable stars ({detectable_count}) is correct, but the specific stars do not match the answer's reasoning. "
                f"Code found: {sorted(detectable_stars_list)}. Expected: {sorted(list(expected_detectable_stars))}.")

    return "Correct"

# Run the check
result = check_correctness_of_answer()
print(result)