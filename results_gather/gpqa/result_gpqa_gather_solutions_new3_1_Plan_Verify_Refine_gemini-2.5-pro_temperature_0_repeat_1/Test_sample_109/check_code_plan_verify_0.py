import math

def check_answer():
    """
    Checks the correctness of the LLM's answer by verifying the observability of each star.
    A star is observable by both ESPRESSO and HIRES if it meets two criteria:
    1. Visibility: Its declination is within the range observable by both observatories.
    2. Brightness: Its apparent V magnitude is brighter than the stricter limit of the two instruments.
    """

    # --- Define Observational Constraints ---

    # 1. Brightness Constraint
    # ESPRESSO limit: V < 17 mag
    # HIRES limit: V < 16 mag
    # The combined, stricter limit for observability by both is V < 16.
    MAGNITUDE_LIMIT = 16.0

    # 2. Visibility (Declination) Constraint
    # Paranal Observatory (ESPRESSO) Latitude: ~ -24.6 deg S
    # Keck Observatory (HIRES) Latitude: ~ +19.8 deg N
    # A star is visible if it rises above the horizon.
    # For Keck (Northern Hemisphere): DEC > Latitude - 90 deg = 19.8 - 90 = -70.2 deg
    # For Paranal (Southern Hemisphere): DEC < 90 deg + Latitude = 90 + (-24.6) = 65.4 deg
    MIN_DEC = -70.2
    MAX_DEC = 65.4

    # --- Define Star Data ---
    stars = {
        "Star1": {
            "dec": -75,
            "m_v": None,
            "M_v": 15.5,
            "dist_pc": 10,
            "E_B_V": 0
        },
        "Star2": {
            "dec": 55,
            "m_v": 16.5,
            "M_v": None,
            "dist_pc": 5,
            "E_B_V": 0
        },
        "Star3": {
            "dec": 48,
            "m_v": 15.5,
            "M_v": None,
            "dist_pc": 15,
            "E_B_V": 0.6
        },
        "Star4": {
            "dec": -48,
            "m_v": None,
            "M_v": 15.5,
            "dist_pc": 10,
            "E_B_V": 0.4
        },
        "Star5": {
            "dec": 60,
            "m_v": None,
            "M_v": 16.5,
            "dist_pc": 5,
            "E_B_V": 0
        }
    }

    # --- Analysis ---
    detectable_stars = []
    analysis_log = []

    for name, data in stars.items():
        # Check 1: Visibility (Declination)
        dec = data["dec"]
        is_visible = MIN_DEC < dec < MAX_DEC
        
        # Check 2: Brightness (Apparent Magnitude)
        m_v = data["m_v"]
        if m_v is None:
            # Calculate apparent magnitude using the distance modulus formula
            # m_v = M_v + 5 * log10(d) - 5 + A_v
            # A_v = 3.1 * E(B-V)
            A_v = 3.1 * data["E_B_V"]
            m_v = data["M_v"] + 5 * math.log10(data["dist_pc"]) - 5 + A_v
        
        is_bright_enough = m_v < MAGNITUDE_LIMIT

        # Log the analysis for each star
        log_entry = (
            f"{name}: DEC={dec}°, m_V={m_v:.2f} mag. "
            f"Visibility check ({MIN_DEC:.1f}° < DEC < {MAX_DEC:.1f}°): {'Pass' if is_visible else 'Fail'}. "
            f"Brightness check (m_V < {MAGNITUDE_LIMIT}): {'Pass' if is_bright_enough else 'Fail'}."
        )
        analysis_log.append(log_entry)

        # Determine if the star is detectable by both
        if is_visible and is_bright_enough:
            detectable_stars.append(name)

    # --- Verify the LLM's Answer ---
    # The provided answer is C, which corresponds to Star3 and Star5.
    llm_answer_choice = "C"
    options = {
        "A": ["Star2", "Star3"],
        "B": ["Star4", "Star5"],
        "C": ["Star3", "Star5"],
        "D": ["Star1", "Star4"]
    }
    
    llm_selected_stars = options.get(llm_answer_choice)

    # Compare the calculated correct stars with the LLM's selection
    if set(detectable_stars) == set(llm_selected_stars):
        return "Correct"
    else:
        reason = "The provided answer is incorrect.\n\n"
        reason += "Here is the step-by-step analysis:\n"
        reason += "1. Visibility Constraint: A star must have a declination between -70.2° and +65.4°.\n"
        reason += "2. Brightness Constraint: A star must have an apparent magnitude V < 16.0.\n\n"
        reason += "Analysis Results for Each Star:\n"
        for log in analysis_log:
            reason += f"- {log}\n"
        reason += "\n"
        reason += f"Based on this analysis, the only stars detectable by both instruments are: {sorted(detectable_stars)}.\n"
        reason += f"The provided answer is '{llm_answer_choice}', which corresponds to the pair {sorted(llm_selected_stars)}.\n"
        reason += f"This is incorrect because the correct pair is {sorted(detectable_stars)}."
        return reason

# Run the check
print(check_answer())