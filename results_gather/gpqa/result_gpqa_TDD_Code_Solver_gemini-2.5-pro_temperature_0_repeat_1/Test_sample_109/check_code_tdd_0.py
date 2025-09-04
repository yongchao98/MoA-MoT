import math

def check_star_observability():
    """
    Checks which stars are observable by both ESPRESSO and HIRES based on
    declination and apparent magnitude constraints.
    """
    # --- Define Observational Constraints ---

    # Stricter magnitude limit (HIRES)
    MAG_LIMIT = 16.0

    # Declination limits based on observatory latitudes
    # Keck Observatory Latitude: ~+19.8 deg
    # Paranal Observatory Latitude: ~-24.6 deg
    # A star at declination 'dec' is visible from latitude 'lat' if lat - 90 <= dec <= lat + 90
    # Keck visible range: [-70.2, 109.8]
    # Paranal visible range: [-114.6, 65.4]
    # The star must be in the intersection of these two ranges.
    DEC_MIN = 19.8 - 90  # -70.2 degrees
    DEC_MAX = -24.6 + 90 # +65.4 degrees

    # --- Star Data ---
    # Note: RA is not needed for this problem.
    stars = {
        "Star1": {"dec": -75, "Mv": 15.5, "dist_pc": 10},
        "Star2": {"dec": 55, "V": 16.5},
        "Star3": {"dec": 48, "V": 15.5},
        "Star4": {"dec": -48, "Mv": 15.5, "dist_pc": 10, "E(B-V)": 0.4},
        "Star5": {"dec": 60, "Mv": 16.5, "dist_pc": 5},
    }

    # --- Analysis ---
    observable_stars = []
    reasons = {}

    for name, data in stars.items():
        # 1. Check Declination
        dec = data["dec"]
        if not (DEC_MIN <= dec <= DEC_MAX):
            reasons[name] = f"Declination {dec}° is outside the observable range [{DEC_MIN:.1f}°, {DEC_MAX:.1f}°]."
            continue

        # 2. Check Apparent Magnitude
        # If apparent magnitude 'V' is given, use it directly.
        if "V" in data:
            v_mag = data["V"]
        # Otherwise, calculate it.
        else:
            mv = data["Mv"]
            dist = data["dist_pc"]
            ebv = data.get("E(B-V)", 0) # Default E(B-V) to 0 if not provided
            
            # Distance modulus: 5 * log10(d) - 5
            distance_modulus = 5 * math.log10(dist) - 5
            # V-band absorption: Av = 3.1 * E(B-V)
            absorption = 3.1 * ebv
            # Apparent magnitude: V = Mv + distance_modulus + Av
            v_mag = mv + distance_modulus + absorption

        if v_mag >= MAG_LIMIT:
            reasons[name] = f"Apparent magnitude V={v_mag:.2f} is not brighter than the {MAG_LIMIT} mag limit."
            continue

        # If both checks pass, the star is observable
        observable_stars.append(name)
        reasons[name] = f"Observable (V={v_mag:.2f}, Dec={dec}°)."

    # --- Verification ---
    # The provided answer 'A' corresponds to "Star3 and Star5"
    expected_answer = {"Star3", "Star5"}
    calculated_answer = set(observable_stars)

    if expected_answer == calculated_answer:
        return "Correct"
    else:
        error_message = "Incorrect. The analysis of the stars is as follows:\n"
        for star, reason in sorted(reasons.items()):
            error_message += f"- {star}: {reason}\n"
        
        error_message += f"\nBased on this, the observable stars are {sorted(list(calculated_answer))}, "
        error_message += f"but the provided answer implies they are {sorted(list(expected_answer))}."
        return error_message

# Run the check and print the result
result = check_star_observability()
print(result)