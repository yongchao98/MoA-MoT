import math

def check_correctness():
    """
    Checks the correctness of the given answer by evaluating each star against the problem's constraints.
    """
    
    # --- Define Observational Constraints ---

    # 1. Visibility (Declination) Constraint
    # Paranal Observatory latitude (approx)
    paranal_lat = -24.6
    # Keck Observatory latitude (approx)
    keck_lat = 19.8
    
    # A star is theoretically visible if its declination (dec) is > (latitude - 90).
    # To be visible from both, it must be in the intersection of their visible ranges.
    # Keck's southern limit: dec > 19.8 - 90 = -70.2 degrees
    # Paranal's northern limit: dec < -24.6 + 90 = 65.4 degrees
    # So, the star's declination must be between -70.2 and +65.4 degrees.
    dec_min_limit = keck_lat - 90
    dec_max_limit = paranal_lat + 90

    # 2. Brightness (Magnitude) Constraint
    # ESPRESSO limit: m_V < 17.0
    # HIRES limit: m_V < 16.0
    # To be detected by BOTH, the star must be brighter than the stricter limit.
    combined_mag_limit = 16.0

    # --- Star Data ---
    # RA is converted from hours to degrees where necessary (1h = 15 deg)
    stars = {
        "Star1": {"DEC": -75, "M_V": 15.5, "d_pc": 10},
        "Star2": {"DEC": 55, "m_V": 16.5},
        "Star3": {"DEC": 48, "m_V": 15.5}, # Apparent magnitude is given, so it's the observed value
        "Star4": {"DEC": -48, "M_V": 15.5, "d_pc": 10, "E(B-V)": 0.4},
        "Star5": {"DEC": 60, "M_V": 16.5, "d_pc": 5},
    }

    # --- Analysis ---
    
    analysis_results = {}
    
    for name, data in stars.items():
        # Step 1: Check declination
        dec = data["DEC"]
        is_dec_observable = dec_min_limit < dec < dec_max_limit
        
        # Step 2: Calculate apparent magnitude (if not given) and check brightness
        apparent_mag = -999 # Default value
        if "m_V" in data:
            apparent_mag = data["m_V"]
        elif "M_V" in data:
            M_V = data["M_V"]
            d_pc = data["d_pc"]
            # Calculate absorption A_V if color excess E(B-V) is provided
            A_V = 3.1 * data.get("E(B-V)", 0)
            # Distance Modulus: m_V = M_V + 5*log10(d) - 5 + A_V
            apparent_mag = M_V + 5 * math.log10(d_pc) - 5 + A_V
            
        is_mag_observable = apparent_mag < combined_mag_limit
        
        # Final conclusion for the star
        analysis_results[name] = is_dec_observable and is_mag_observable

    # --- Verify the Answer ---
    
    given_answer_option = "D"
    options = {
        "A": ["Star2", "Star3"],
        "B": ["Star4", "Star5"],
        "C": ["Star1", "Star4"],
        "D": ["Star3", "Star5"],
    }
    
    # Find the set of stars that are actually observable
    correctly_observable_stars = {name for name, is_observable in analysis_results.items() if is_observable}
    
    # Get the set of stars from the given answer
    stars_in_given_answer = set(options[given_answer_option])
    
    if correctly_observable_stars == stars_in_given_answer:
        return "Correct"
    else:
        # Generate a detailed reason for the incorrectness
        reason = "The provided answer is incorrect.\n\n"
        reason += "Here is the detailed analysis:\n"
        reason += f"Observability requires: Declination between {dec_min_limit:.1f}° and {dec_max_limit:.1f}°, AND Apparent Magnitude < {combined_mag_limit} mag.\n\n"
        
        # Star 1
        dec1 = stars["Star1"]["DEC"]
        m1 = stars["Star1"]["M_V"] + 5*math.log10(stars["Star1"]["d_pc"]) - 5
        reason += f"Star1: DEC = {dec1}°. Fails declination check (not > {dec_min_limit:.1f}°). Apparent mag = {m1:.2f}. Observable: {analysis_results['Star1']}.\n"
        
        # Star 2
        dec2 = stars["Star2"]["DEC"]
        m2 = stars["Star2"]["m_V"]
        reason += f"Star2: DEC = {dec2}°. OK. Apparent mag = {m2:.2f}. Fails magnitude check (not < {combined_mag_limit}). Observable: {analysis_results['Star2']}.\n"
        
        # Star 3
        dec3 = stars["Star3"]["DEC"]
        m3 = stars["Star3"]["m_V"]
        reason += f"Star3: DEC = {dec3}°. OK. Apparent mag = {m3:.2f}. OK. Observable: {analysis_results['Star3']}.\n"
        
        # Star 4
        dec4 = stars["Star4"]["DEC"]
        m4 = stars["Star4"]["M_V"] + 5*math.log10(stars["Star4"]["d_pc"]) - 5 + 3.1*stars["Star4"]["E(B-V)"]
        reason += f"Star4: DEC = {dec4}°. OK. Apparent mag = {m4:.2f}. Fails magnitude check (not < {combined_mag_limit}). Observable: {analysis_results['Star4']}.\n"
        
        # Star 5
        dec5 = stars["Star5"]["DEC"]
        m5 = stars["Star5"]["M_V"] + 5*math.log10(stars["Star5"]["d_pc"]) - 5
        reason += f"Star5: DEC = {dec5}°. OK. Apparent mag = {m5:.2f}. OK. Observable: {analysis_results['Star5']}.\n"
        
        reason += f"\nConclusion: The only stars observable by both are {sorted(list(correctly_observable_stars))}.\n"
        reason += f"The given answer '{given_answer_option}' ({sorted(list(stars_in_given_answer))}) is not the correct set."
        
        return reason

# Run the check and print the result
result = check_correctness()
print(result)