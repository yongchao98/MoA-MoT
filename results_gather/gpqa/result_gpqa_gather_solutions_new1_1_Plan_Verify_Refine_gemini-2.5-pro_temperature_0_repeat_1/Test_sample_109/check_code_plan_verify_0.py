import math

def check_answer():
    """
    Checks which stars are observable by both ESPRESSO and HIRES based on brightness and visibility constraints.
    """

    # --- Define Constraints ---
    # Brightness: Must be brighter than the stricter limit (HIRES)
    MAGNITUDE_LIMIT = 16.0

    # Visibility: Must be in the declination range visible from both hemispheres
    # Keck Observatory (Latitude ~+19.8 N) can see DEC > 19.8 - 90 = -70.2 deg
    # Paranal Observatory (Latitude ~-24.6 S) can see DEC < -24.6 + 90 = +65.4 deg
    DEC_LOWER_LIMIT = -70.2
    DEC_UPPER_LIMIT = 65.4

    # --- Star Data ---
    # RA is not needed for this analysis.
    stars = [
        {'name': 'Star1', 'dec': -75, 'M_v': 15.5, 'd_pc': 10},
        {'name': 'Star2', 'dec': 55, 'V_app': 16.5, 'd_pc': 5},
        {'name': 'Star3', 'dec': 48, 'V_app': 15.5, 'E_BV': 0.6, 'd_pc': 15},
        {'name': 'Star4', 'dec': -48, 'M_v': 15.5, 'd_pc': 10, 'E_BV': 0.4},
        {'name': 'Star5', 'dec': 60, 'M_v': 16.5, 'd_pc': 5},
    ]

    # --- Helper function to calculate apparent magnitude ---
    def get_apparent_magnitude(star):
        # If apparent magnitude is given directly, use it.
        if 'V_app' in star:
            return star['V_app']
        
        # Otherwise, calculate it from absolute magnitude, distance, and extinction.
        M_v = star['M_v']
        d_pc = star['d_pc']
        A_v = 0
        if 'E_BV' in star:
            A_v = 3.1 * star['E_BV']
        
        # Distance modulus formula: V = M_v + 5 * log10(d_pc / 10) + A_v
        V_app = M_v + 5 * math.log10(d_pc / 10) + A_v
        return V_app

    # --- Analysis ---
    observable_stars = []
    analysis_log = "--- Star-by-Star Analysis ---\n"
    for star in stars:
        name = star['name']
        dec = star['dec']
        v_mag = get_apparent_magnitude(star)

        # Check conditions
        mag_ok = v_mag < MAGNITUDE_LIMIT
        dec_ok = DEC_LOWER_LIMIT < dec < DEC_UPPER_LIMIT

        analysis_log += f"{name}:\n"
        analysis_log += f"  - Apparent Magnitude (V): {v_mag:.2f}. Condition V < {MAGNITUDE_LIMIT}: {'Pass' if mag_ok else 'Fail'}\n"
        analysis_log += f"  - Declination (DEC): {dec}°. Condition {DEC_LOWER_LIMIT:.1f}° < DEC < {DEC_UPPER_LIMIT:.1f}°: {'Pass' if dec_ok else 'Fail'}\n"
        
        if mag_ok and dec_ok:
            observable_stars.append(name)
            analysis_log += "  - Result: Observable by both.\n"
        else:
            analysis_log += "  - Result: Not observable by both.\n"

    # --- Verification ---
    # The final answer from the LLM is B, which corresponds to Star3 and Star5.
    expected_stars = {'Star3', 'Star5'}
    
    if set(observable_stars) == expected_stars:
        return "Correct"
    else:
        reason = "The provided answer is incorrect.\n\n"
        reason += f"My analysis concluded that the observable stars are {sorted(observable_stars)}.\n"
        reason += f"The provided answer corresponds to {sorted(list(expected_stars))}.\n\n"
        reason += analysis_log
        return reason

# Run the check
result = check_answer()
print(result)