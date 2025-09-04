import math

def check_star_observability():
    """
    This function checks which stars from a predefined list can be observed by
    both ESPRESSO (Paranal) and HIRES (Keck) spectrographs based on magnitude
    and location constraints. It then compares the result to the expected answer.
    """

    # --- Define Constraints ---

    # Paranal Observatory (ESPRESSO)
    lat_paranal = -24.6  # Latitude in degrees
    mag_limit_espresso = 17.0  # Apparent V magnitude limit

    # W. M. Keck Observatory (HIRES)
    lat_keck = 19.8  # Latitude in degrees
    mag_limit_hires = 16.0  # Apparent V magnitude limit

    # --- Define Star Data ---
    # Note: RA is not needed for the calculation as per the problem statement.
    stars = {
        "Star1": {
            "dec": -75.0, "Mv": 15.5, "d_pc": 10.0, "mv": None, "E(B-V)": 0.0
        },
        "Star2": {
            "dec": 55.0, "mv": 16.5, "Mv": None, "d_pc": 5.0, "E(B-V)": 0.0
        },
        "Star3": {
            "dec": 48.0, "mv": 15.5, "Mv": None, "d_pc": 15.0, "E(B-V)": 0.6
        },
        "Star4": {
            "dec": -48.0, "Mv": 15.5, "d_pc": 10.0, "mv": None, "E(B-V)": 0.4
        },
        "Star5": {
            "dec": 60.0, "Mv": 16.5, "d_pc": 5.0, "mv": None, "E(B-V)": 0.0
        }
    }

    # --- Calculation and Verification ---

    # Visibility condition: A star is visible if it rises above the horizon.
    # For Keck (Northern Hemisphere): DEC > Lat - 90
    # For Paranal (Southern Hemisphere): DEC < Lat + 90
    keck_dec_min = lat_keck - 90.0
    paranal_dec_max = lat_paranal + 90.0

    observable_by_both = []
    analysis_log = {}

    for name, data in stars.items():
        # Step 1: Determine the apparent V magnitude (mv)
        # If mv is given directly, use it. Otherwise, calculate it.
        current_mv = data["mv"]
        if current_mv is None:
            # Formula: mv = Mv - 5 + 5*log10(d) + Av
            # where Av = 3.1 * E(B-V)
            distance_modulus = 5 * math.log10(data["d_pc"]) - 5
            extinction = 3.1 * data["E(B-V)"]
            current_mv = data["Mv"] + distance_modulus + extinction
        
        # Step 2: Check magnitude constraints for both observatories
        is_bright_enough_for_espresso = current_mv < mag_limit_espresso
        is_bright_enough_for_hires = current_mv < mag_limit_hires

        if not is_bright_enough_for_hires:
            analysis_log[name] = f"Fails HIRES magnitude limit. Calculated mv ({current_mv:.2f}) is not < {mag_limit_hires}."
            continue
        if not is_bright_enough_for_espresso:
            analysis_log[name] = f"Fails ESPRESSO magnitude limit. Calculated mv ({current_mv:.2f}) is not < {mag_limit_espresso}."
            continue

        # Step 3: Check visibility (declination) constraints for both observatories
        dec = data["dec"]
        is_visible_from_keck = dec > keck_dec_min
        is_visible_from_paranal = dec < paranal_dec_max

        if not is_visible_from_keck:
            analysis_log[name] = f"Not visible from Keck. Declination ({dec}°) is not > {keck_dec_min:.1f}°."
            continue
        if not is_visible_from_paranal:
            analysis_log[name] = f"Not visible from Paranal. Declination ({dec}°) is not < {paranal_dec_max:.1f}°."
            continue
        
        # If all checks pass, the star is observable by both
        observable_by_both.append(name)
        analysis_log[name] = "Passes all checks."

    # Step 4: Compare the result with the LLM's answer (B: Star3 and Star5)
    expected_stars = {"Star3", "Star5"}
    
    if set(observable_by_both) == expected_stars:
        return "Correct"
    else:
        error_message = "Incorrect. The final list of observable stars does not match the expected answer.\n"
        error_message += f"Expected: {sorted(list(expected_stars))}\n"
        error_message += f"Calculated: {sorted(observable_by_both)}\n\n"
        error_message += "Detailed analysis per star:\n"
        for name, result in sorted(analysis_log.items()):
            error_message += f"- {name}: {result}\n"
        return error_message.strip()

# Run the check and print the result
result = check_star_observability()
print(result)