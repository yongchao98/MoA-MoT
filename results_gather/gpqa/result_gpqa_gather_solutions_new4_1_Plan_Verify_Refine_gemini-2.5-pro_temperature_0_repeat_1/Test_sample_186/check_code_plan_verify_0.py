import math

def check_correctness_of_answer():
    """
    This function checks the correctness of the answer to the ESPRESSO detectability question.
    It verifies the two main criteria for each star:
    1. Visibility from Paranal Observatory.
    2. Brightness to achieve S/N >= 10 in a 1-hour exposure.
    """

    # --- Step 1: Define constants and detection criteria ---

    # Paranal Observatory Latitude in degrees
    PARANAL_LATITUDE = -24.6

    # Visibility limit: A star is not visible if its declination is too far north.
    # The maximum declination visible is approximately 90 - |latitude|.
    VISIBILITY_DEC_LIMIT = 90 - abs(PARANAL_LATITUDE)

    # Brightness limit: Apparent magnitude (V) must be <= 17.0 for S/N >= 10 in 1hr.
    # This is based on the ESO ESPRESSO performance table for a 3600s exposure.
    LIMITING_MAGNITUDE = 17.0

    # --- Step 2: Define the list of stars to be evaluated ---

    stars = [
        {
            "name": "a) Canopus",
            "apparent_mag": -0.74,
            "declination": -52.7,
            "absolute_mag": None,
            "distance_pc": None,
        },
        {
            "name": "b) Polaris",
            "apparent_mag": 1.98,
            "declination": 89.25,
            "absolute_mag": None,
            "distance_pc": None,
        },
        {
            "name": "c) Star at 10 pc",
            "apparent_mag": None, # To be calculated
            "declination": 0.0,
            "absolute_mag": 15.0,
            "distance_pc": 10.0,
        },
        {
            "name": "d) Star at 200 pc",
            "apparent_mag": None,
            "declination": 0.0,
            "absolute_mag": 15.0,
            "distance_pc": 200.0,
        },
        {
            "name": "e) Star at 5 pc",
            "apparent_mag": None,
            "declination": 0.0,
            "absolute_mag": 15.0,
            "distance_pc": 5.0,
        },
        {
            "name": "f) Star at 50 pc",
            "apparent_mag": None,
            "declination": 0.0,
            "absolute_mag": 15.0,
            "distance_pc": 50.0,
        },
    ]

    # --- Step 3: Helper function to calculate apparent magnitude ---

    def calculate_apparent_magnitude(Mv, d):
        """Calculates apparent magnitude (V) from absolute magnitude (Mv) and distance (d) in parsecs."""
        if d <= 0:
            return float('inf')  # Invalid distance
        # Formula: V = Mv - 5 + 5 * log10(d)
        return Mv - 5 + 5 * math.log10(d)

    # --- Step 4: Evaluate each star and count detectable ones ---

    detectable_count = 0
    evaluation_log = []

    for star in stars:
        # Condition 1: Check visibility
        is_visible = star["declination"] < VISIBILITY_DEC_LIMIT
        if not is_visible:
            evaluation_log.append(
                f"Star '{star['name']}' is NOT detectable. Reason: Not visible from Paranal (Declination {star['declination']:.2f}° > {VISIBILITY_DEC_LIMIT:.2f}°)."
            )
            continue

        # Get or calculate apparent magnitude
        if star["apparent_mag"] is not None:
            apparent_mag = star["apparent_mag"]
        else:
            apparent_mag = calculate_apparent_magnitude(star["absolute_mag"], star["distance_pc"])

        # Condition 2: Check brightness
        is_bright_enough = apparent_mag <= LIMITING_MAGNITUDE
        if not is_bright_enough:
            evaluation_log.append(
                f"Star '{star['name']}' is NOT detectable. Reason: Too faint (Apparent Mag {apparent_mag:.2f} > {LIMITING_MAGNITUDE})."
            )
            continue

        # If both conditions are met, the star is detectable
        detectable_count += 1
        evaluation_log.append(
            f"Star '{star['name']}' IS detectable. (Visible and Apparent Mag {apparent_mag:.2f} <= {LIMITING_MAGNITUDE})."
        )

    # --- Step 5: Compare with the provided answer ---
    # The provided answer is D, which corresponds to 3 detectable stars.
    expected_count = 3

    if detectable_count == expected_count:
        return "Correct"
    else:
        error_message = (
            f"Incorrect. The code's analysis found {detectable_count} detectable stars, "
            f"but the provided answer implies {expected_count}.\n\n"
            "Detailed Evaluation Log:\n"
        )
        for log_entry in evaluation_log:
            error_message += f"- {log_entry}\n"
        return error_message

# Execute the check and print the result.
print(check_correctness_of_answer())