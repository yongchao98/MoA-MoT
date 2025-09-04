import math

def check_answer():
    """
    This function checks the correctness of the provided answer by verifying the two main constraints:
    1. Visibility: The star must be observable from both observatories.
    2. Brightness: The star's apparent magnitude must be brighter than the stricter limit of the two instruments.
    """

    # --- Define Constants and Constraints ---

    # Latitudes of the observatories
    LAT_PARANAL = -24.6  # ESPRESSO, Southern Hemisphere
    LAT_KECK = 19.8      # HIRES, Northern Hemisphere

    # Magnitude limits for the spectrographs
    MAG_LIMIT_ESPRESSO = 17.0
    MAG_LIMIT_HIRES = 16.0

    # Derived constraints based on the problem statement
    # For a star to be visible from both, its declination must be in the overlapping range.
    # Keck can see south to DEC > (LAT_KECK - 90)
    # Paranal can see north to DEC < (LAT_PARANAL + 90)
    DEC_MIN = LAT_KECK - 90.0  # Approx -70.2 degrees
    DEC_MAX = LAT_PARANAL + 90.0 # Approx +65.4 degrees

    # For a star to be detected by both, it must be brighter than the stricter limit.
    MAG_LIMIT_COMBINED = min(MAG_LIMIT_ESPRESSO, MAG_LIMIT_HIRES)

    # --- Star Data ---
    # RA is not needed for the calculation.
    # If m_V (apparent magnitude) is given, it's the final observed value.
    # If M_V (absolute magnitude) is given, m_V must be calculated.
    stars = [
        {'name': 'Star1', 'DEC': -75, 'M_V': 15.5, 'd': 10},
        {'name': 'Star2', 'DEC': 55, 'm_V': 16.5, 'd': 5},
        {'name': 'Star3', 'DEC': 48, 'm_V': 15.5, 'E(B-V)': 0.6, 'd': 15},
        {'name': 'Star4', 'DEC': -48, 'M_V': 15.5, 'E(B-V)': 0.4, 'd': 10},
        {'name': 'Star5', 'DEC': 60, 'M_V': 16.5, 'd': 5},
    ]

    # --- Analysis ---
    observable_stars = []
    reasons = {}

    for star in stars:
        name = star['name']
        
        # 1. Check Visibility (Declination)
        dec = star['DEC']
        is_visible = DEC_MIN < dec < DEC_MAX
        
        if not is_visible:
            reasons[name] = f"Fails visibility constraint. Its declination {dec}° is outside the observable range ({DEC_MIN:.1f}° to {DEC_MAX:.1f}°)."
            continue

        # 2. Check Brightness (Apparent Magnitude)
        # If apparent magnitude 'm_V' is given, use it directly.
        if 'm_V' in star:
            m_v = star['m_V']
        else:
            # Otherwise, calculate it from absolute magnitude 'M_V'.
            # Formula: m_V = M_V + 5 * log10(d/10) + A_V
            M_V = star['M_V']
            d = star['d']
            
            # Calculate extinction A_V if E(B-V) is provided
            A_V = 3.1 * star.get('E(B-V)', 0)
            
            m_v = M_V + 5 * math.log10(d / 10) + A_V

        is_bright_enough = m_v < MAG_LIMIT_COMBINED

        if not is_bright_enough:
            reasons[name] = f"Fails brightness constraint. Its apparent magnitude is {m_v:.2f}, which is not brighter than the required {MAG_LIMIT_COMBINED}."
            continue
            
        # If both checks pass, the star is observable by both.
        observable_stars.append(name)
        reasons[name] = f"Passes both constraints (DEC={dec}°, V_mag={m_v:.2f})."

    # --- Final Verification ---
    # The provided answer is 'A', which corresponds to 'Star3 and Star5'.
    expected_stars = {'Star3', 'Star5'}
    
    # Convert the list of found stars to a set for easy comparison.
    found_stars_set = set(observable_stars)

    if found_stars_set == expected_stars:
        return "Correct"
    else:
        error_message = "Incorrect. The analysis of which stars are observable is wrong.\n"
        error_message += f"Expected observable stars: {sorted(list(expected_stars))}\n"
        error_message += f"Calculated observable stars: {sorted(list(found_stars_set))}\n\n"
        error_message += "Detailed breakdown:\n"
        for star in stars:
            name = star['name']
            error_message += f"- {name}: {reasons[name]}\n"
        
        return error_message

# Run the check
result = check_answer()
print(result)