import math

def check_correctness():
    """
    This function verifies the solution to the star detection problem.
    It calculates the number of detectable stars based on the problem's constraints
    and compares it to the provided answer.
    """

    # --- Define Constants and Formulas ---
    PARANAL_LATITUDE = -24.6275  # Latitude of Paranal Observatory in degrees

    # Calculate the limiting magnitude for ESPRESSO in 1-UT mode
    # S/N=10 for V=20 in 1hr in 4-UT mode. 1-UT needs 4x the flux.
    # Δm = 2.5 * log10(flux_ratio)
    LIMITING_MAGNITUDE = 20 - 2.5 * math.log10(4)

    # --- Helper Functions ---
    def calculate_apparent_magnitude(M_V, d_pc):
        """Calculates apparent magnitude (m) from absolute magnitude (M) and distance (d)."""
        if d_pc <= 0:
            return float('inf')
        # Distance modulus formula: m = M + 5*log10(d) - 5
        return M_V + 5 * math.log10(d_pc) - 5

    def is_observable(dec_deg, lat_deg):
        """
        Checks if a star is ever visible from a given latitude.
        A star is observable if it is not circumpolar around the opposite celestial pole.
        For a southern observatory, this means dec < (90 - abs(lat)).
        """
        return dec_deg < (90 - abs(lat_deg))

    # --- Star Data from the Question ---
    stars = [
        {'name': 'a) Canopus', 'm_V': -0.74, 'dec': -52.7},
        {'name': 'b) Polaris', 'm_V': 1.98, 'dec': 89.26},
        {'name': 'c) Star 10pc', 'M_V': 15, 'dist_pc': 10, 'dec': 0},
        {'name': 'd) Star 200pc', 'M_V': 15, 'dist_pc': 200, 'dec': 0},
        {'name': 'e) Star 5pc', 'M_V': 15, 'dist_pc': 5, 'dec': 0},
        {'name': 'f) Star 50pc', 'M_V': 15, 'dist_pc': 50, 'dec': 0},
    ]

    # --- Main Logic: Count Detectable Stars ---
    detectable_stars_count = 0
    analysis_log = []

    for star in stars:
        # Determine apparent magnitude
        if 'm_V' in star:
            m_V = star['m_V']
        else:
            m_V = calculate_apparent_magnitude(star['M_V'], star['dist_pc'])

        # Check conditions
        observable = is_observable(star['dec'], PARANAL_LATITUDE)
        bright_enough = (m_V <= LIMITING_MAGNITUDE)
        
        is_detectable = observable and bright_enough
        if is_detectable:
            detectable_stars_count += 1
        
        analysis_log.append(
            f"Star: {star['name']}\n"
            f"  - Apparent Mag (V): {m_V:.3f}\n"
            f"  - Declination (dec): {star['dec']:.2f}°\n"
            f"  - Is Observable? {observable} (dec < {90 - abs(PARANAL_LATITUDE):.2f}°)\n"
            f"  - Is Bright Enough? {bright_enough} (V <= {LIMITING_MAGNITUDE:.3f})\n"
            f"  - DETECTABLE: {is_detectable}\n"
        )

    # --- Verification ---
    # The provided answer is B, which corresponds to a count of 4.
    expected_count = 4
    
    if detectable_stars_count == expected_count:
        return "Correct"
    else:
        error_message = (
            f"The answer is incorrect. The calculated number of detectable stars is {detectable_stars_count}, "
            f"but the provided answer 'B' implies {expected_count}.\n\n"
            "Here is the detailed analysis:\n"
            "-------------------------------------\n"
        )
        error_message += "\n".join(analysis_log)
        return error_message

# Run the check and print the result
result = check_correctness()
print(result)