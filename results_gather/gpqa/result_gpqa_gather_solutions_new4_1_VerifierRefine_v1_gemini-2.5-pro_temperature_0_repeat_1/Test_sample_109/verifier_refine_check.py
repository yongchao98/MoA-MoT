import math

def check_answer():
    """
    This function checks the correctness of the provided answer by re-evaluating the conditions for each star.
    """
    
    # --- Define Observational Constraints ---
    
    # Brightness limits
    mag_limit_espresso = 17.0
    mag_limit_hires = 16.0
    # The combined limit is the stricter of the two
    combined_mag_limit = min(mag_limit_espresso, mag_limit_hires)

    # Positional (Declination) limits based on observatory latitudes
    # Paranal Observatory (ESPRESSO): Latitude ~ -24.6° S
    # W. M. Keck Observatory (HIRES): Latitude ~ +19.8° N
    # The question asks to disregard pointing limits, so we check if the star ever rises above the horizon.
    # For Keck (Northern Hemisphere), the southern limit is DEC > Lat - 90°
    dec_limit_south = 19.8 - 90.0  # approx -70.2°
    # For Paranal (Southern Hemisphere), the northern limit is DEC < Lat + 90°
    dec_limit_north = -24.6 + 90.0 # approx +65.4°

    # --- Define Star Data ---
    
    stars = {
        "Star1": {"RA": 15, "DEC": -75, "M_V": 15.5, "d": 10, "E(B-V)": 0},
        "Star2": {"RA": 30, "DEC": 55, "m_V": 16.5, "d": 5},
        "Star3": {"RA_h": 11, "DEC": 48, "m_V": 15.5, "E(B-V)": 0.6, "d": 15},
        "Star4": {"RA": 85, "DEC": -48, "M_V": 15.5, "E(B-V)": 0.4, "d": 10},
        "Star5": {"RA_h": 10, "DEC": 60, "M_V": 16.5, "d": 5, "E(B-V)": 0}
    }

    observable_stars = []
    reasons = []

    # --- Analyze Each Star ---
    
    for name, data in stars.items():
        # 1. Check Positional Visibility (Declination)
        dec = data["DEC"]
        is_visible_positionally = dec_limit_south < dec < dec_limit_north
        
        if not is_visible_positionally:
            reasons.append(f"{name}: FAILED visibility check. Declination {dec}° is outside the common observable range ({dec_limit_south:.1f}° to {dec_limit_north:.1f}°).")
            continue

        # 2. Check Brightness (Apparent Magnitude)
        apparent_mag = 0
        # If apparent magnitude is given directly, use it.
        if "m_V" in data:
            apparent_mag = data["m_V"]
            # A careful point: For Star3, the problem gives an apparent magnitude. This is the final observed brightness,
            # so the E(B-V) and distance are extraneous for the magnitude check.
        else:
            # Otherwise, calculate it using the distance modulus formula: m_V = M_V + 5*log10(d) - 5 + A_V
            M_V = data["M_V"]
            d = data["d"]
            # Calculate extinction A_V = 3.1 * E(B-V)
            A_V = 3.1 * data.get("E(B-V)", 0)
            apparent_mag = M_V + 5 * math.log10(d) - 5 + A_V

        is_bright_enough = apparent_mag < combined_mag_limit

        if not is_bright_enough:
            reasons.append(f"{name}: FAILED brightness check. Apparent magnitude {apparent_mag:.2f} is not brighter than the required limit of {combined_mag_limit}.")
            continue
            
        # If both checks pass, the star is observable
        observable_stars.append(name)
        reasons.append(f"{name}: PASSED. Declination {dec}° is visible. Apparent magnitude {apparent_mag:.2f} is brighter than {combined_mag_limit}.")

    # --- Final Verification ---
    
    # The provided answer is D, which corresponds to Star3 and Star5.
    expected_answer = ["Star3", "Star5"]
    
    # Sort both lists to ensure comparison is order-independent
    observable_stars.sort()
    expected_answer.sort()

    if observable_stars == expected_answer:
        return "Correct"
    else:
        error_message = "Incorrect. The analysis leads to a different set of observable stars.\n"
        error_message += "Analysis details:\n"
        for reason in reasons:
            error_message += f"- {reason}\n"
        error_message += f"\nCalculated observable stars: {observable_stars}\n"
        error_message += f"Expected observable stars (from answer D): {expected_answer}"
        return error_message

# Run the check
result = check_answer()
print(result)