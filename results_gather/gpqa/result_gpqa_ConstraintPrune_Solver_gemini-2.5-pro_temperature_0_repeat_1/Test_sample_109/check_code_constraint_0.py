import math

def check_answer():
    """
    This function checks the correctness of the provided answer by verifying
    which stars can be observed by both ESPRESSO and HIRES spectrographs.

    The verification is based on two main constraints:
    1. Visibility: The star's declination must be within the observable range of both observatories.
    2. Brightness: The star's apparent V magnitude must be brighter than the stricter limit of the two instruments.
    """

    # --- Define Observational Constraints ---

    # 1. Visibility (Declination) Constraint
    # Latitude of Paranal Observatory (ESPRESSO): ~ -24.6 degrees
    # Latitude of Keck Observatory (HIRES): ~ +19.8 degrees
    # A star is visible if its declination (DEC) is between:
    # Keck's southern limit: 19.8 - 90 = -70.2 degrees
    # Paranal's northern limit: 90 - 24.6 = 65.4 degrees
    MIN_DEC = -70.2
    MAX_DEC = 65.4

    # 2. Brightness (Apparent Magnitude) Constraint
    # ESPRESSO limit: V < 17 mag
    # HIRES limit: V < 16 mag
    # To be observed by BOTH, the star must satisfy the stricter limit.
    MAX_V_MAG = 16.0

    # --- Star Data ---
    # A dictionary to hold the properties of each star.
    # Note: RA is not needed for the calculation but is included for completeness.
    stars = {
        "Star1": {"DEC": -75, "M_V": 15.5, "dist_pc": 10, "E(B-V)": 0.0},
        "Star2": {"DEC": 55, "V": 16.5},
        "Star3": {"DEC": 48, "V": 15.5},
        "Star4": {"DEC": -48, "M_V": 15.5, "dist_pc": 10, "E(B-V)": 0.4},
        "Star5": {"DEC": 60, "M_V": 16.5, "dist_pc": 5, "E(B-V)": 0.0},
    }

    # --- Analysis ---
    observable_stars = []
    analysis_log = []

    for name, data in stars.items():
        # Step 1: Check Declination Constraint
        dec = data["DEC"]
        dec_pass = MIN_DEC < dec < MAX_DEC

        if not dec_pass:
            analysis_log.append(f"{name}: FAILED. Declination {dec} is outside the observable range ({MIN_DEC} to {MAX_DEC}).")
            continue

        # Step 2: Check Apparent Magnitude Constraint
        # Calculate apparent magnitude 'V' if not provided directly.
        if "V" in data:
            v_mag = data["V"]
        else:
            # Use the distance modulus formula: V = M_V + 5*log10(d) - 5 + A_V
            # where A_V = 3.1 * E(B-V)
            m_v = data["M_V"]
            dist_pc = data["dist_pc"]
            e_b_v = data["E(B-V)"]
            a_v = 3.1 * e_b_v
            v_mag = m_v + 5 * math.log10(dist_pc) - 5 + a_v

        mag_pass = v_mag < MAX_V_MAG

        if mag_pass:
            analysis_log.append(f"{name}: PASSED. Declination ({dec}) is valid and Apparent Magnitude ({v_mag:.2f}) is < {MAX_V_MAG}.")
            observable_stars.append(name)
        else:
            analysis_log.append(f"{name}: FAILED. Apparent Magnitude ({v_mag:.2f}) is not < {MAX_V_MAG}.")

    # --- Final Verification ---
    # The provided answer 'A' corresponds to "Star3 and Star5".
    expected_stars = ["Star3", "Star5"]

    # Sort lists to ensure comparison is order-independent
    if sorted(observable_stars) == sorted(expected_stars):
        return "Correct"
    else:
        error_message = "Incorrect. The calculated observable stars do not match the provided answer.\n"
        error_message += "Provided answer implies Star3 and Star5 are observable.\n"
        error_message += "My analysis found the following stars to be observable: "
        error_message += f"{observable_stars if observable_stars else 'None'}\n\n"
        error_message += "Detailed Analysis Log:\n"
        error_message += "\n".join(analysis_log)
        return error_message

# Execute the check and print the result
result = check_answer()
print(result)