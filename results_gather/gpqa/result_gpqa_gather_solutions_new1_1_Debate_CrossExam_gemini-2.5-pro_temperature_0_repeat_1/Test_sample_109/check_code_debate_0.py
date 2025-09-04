import math

def check_answer():
    """
    Checks the correctness of the LLM's answer by verifying the observability
    of each star based on location and brightness constraints.
    """

    # --- Define Observational Constraints ---
    # Brightness: Must be brighter than the stricter limit (HIRES)
    MAGNITUDE_LIMIT = 16.0

    # Location: Must be in the overlapping declination range
    # Paranal (Lat ~ -24.6 S) can see DEC < -24.6 + 90 = +65.4
    # Keck (Lat ~ +19.8 N) can see DEC > 19.8 - 90 = -70.2
    DEC_UPPER_LIMIT = 65.4
    DEC_LOWER_LIMIT = -70.2

    # --- Star Data ---
    stars = {
        "Star1": {"dec": -75, "Mv": 15.5, "d_pc": 10, "E_BV": 0},
        "Star2": {"dec": 55, "V": 16.5},
        "Star3": {"dec": 48, "V": 15.5}, # Apparent V mag is given directly
        "Star4": {"dec": -48, "Mv": 15.5, "d_pc": 10, "E_BV": 0.4},
        "Star5": {"dec": 60, "Mv": 16.5, "d_pc": 5, "E_BV": 0},
    }

    observable_stars = []
    reasons = {}

    # --- Analysis Loop ---
    for name, data in stars.items():
        # 1. Location Check
        dec = data["dec"]
        is_location_ok = DEC_LOWER_LIMIT < dec < DEC_UPPER_LIMIT
        if not is_location_ok:
            reasons[name] = f"Fails location check. DEC={dec} is outside the observable range ({DEC_LOWER_LIMIT:.1f} to {DEC_UPPER_LIMIT:.1f})."
            continue

        # 2. Brightness Check
        apparent_mag = 0
        if "V" in data:
            apparent_mag = data["V"]
        else:
            # Calculate apparent magnitude from absolute magnitude and distance
            Mv = data["Mv"]
            d_pc = data["d_pc"]
            E_BV = data["E_BV"]
            
            # Calculate absorption Av = 3.1 * E(B-V)
            Av = 3.1 * E_BV
            
            # Distance Modulus: V = Mv + 5*log10(d/10) + Av
            # Note: log10(d/10) can be written as log10(d) - 1
            distance_modulus = 5 * (math.log10(d_pc) - 1)
            apparent_mag = Mv + distance_modulus + Av

        is_brightness_ok = apparent_mag < MAGNITUDE_LIMIT
        if not is_brightness_ok:
            reasons[name] = f"Fails brightness check. Apparent V mag = {apparent_mag:.2f}, which is not brighter than {MAGNITUDE_LIMIT}."
            continue
            
        # If both checks pass, the star is observable
        observable_stars.append(name)
        reasons[name] = f"Passes both checks. Apparent V mag = {apparent_mag:.2f}, DEC = {dec}."

    # --- Final Verification ---
    # The provided answer is 'A', which corresponds to Star3 and Star5.
    expected_answer = ["Star3", "Star5"]
    
    # Sort for consistent comparison
    observable_stars.sort()
    expected_answer.sort()

    if observable_stars == expected_answer:
        return "Correct"
    else:
        error_message = "Incorrect. The analysis shows a different set of observable stars.\n"
        error_message += f"Expected observable stars: {expected_answer}\n"
        error_message += f"Calculated observable stars: {observable_stars}\n\n"
        error_message += "Detailed breakdown:\n"
        for star, reason in sorted(reasons.items()):
            error_message += f"- {star}: {reason}\n"
        return error_message

# Run the check and print the result
result = check_answer()
print(result)