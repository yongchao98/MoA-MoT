import math

def check_star_observability():
    """
    Checks which stars can be observed by both ESPRESSO and HIRES based on
    visibility and brightness constraints.
    """
    # --- Define Constraints ---
    # Brightness limit (stricter of the two)
    MAGNITUDE_LIMIT = 16.0

    # Visibility limits based on observatory latitudes
    # Paranal Lat ~ -24.6 deg -> Northern DEC limit = -24.6 + 90 = 65.4
    # Keck Lat ~ 19.8 deg -> Southern DEC limit = 19.8 - 90 = -70.2
    DEC_MAX = 65.4
    DEC_MIN = -70.2

    # --- Star Data ---
    stars = [
        {'name': 'Star1', 'dec': -75, 'M_v': 15.5, 'd': 10, 'E_BV': 0, 'm_v': None},
        {'name': 'Star2', 'dec': 55, 'm_v': 16.5},
        {'name': 'Star3', 'dec': 48, 'm_v': 15.5}, # Apparent magnitude is given, other data is extraneous
        {'name': 'Star4', 'dec': -48, 'M_v': 15.5, 'd': 10, 'E_BV': 0.4, 'm_v': None},
        {'name': 'Star5', 'dec': 60, 'M_v': 16.5, 'd': 5, 'E_BV': 0, 'm_v': None},
    ]

    observable_stars = []
    analysis_log = []

    for star in stars:
        name = star['name']
        dec = star['dec']
        
        # 1. Check Visibility
        is_visible = DEC_MIN < dec < DEC_MAX
        
        if not is_visible:
            analysis_log.append(f"{name}: FAILED visibility check. DEC={dec} is outside the range ({DEC_MIN}, {DEC_MAX}).")
            continue

        # 2. Check Brightness
        apparent_magnitude = 0
        if star.get('m_v') is not None:
            apparent_magnitude = star['m_v']
        else:
            # Calculate apparent magnitude
            M_v = star['M_v']
            d = star['d']
            E_BV = star['E_BV']
            A_v = 3.1 * E_BV
            # m_v = M_v + 5 * log10(d/10) + A_v
            apparent_magnitude = M_v + 5 * math.log10(d / 10) + A_v

        is_bright_enough = apparent_magnitude < MAGNITUDE_LIMIT

        if not is_bright_enough:
            analysis_log.append(f"{name}: FAILED brightness check. Apparent magnitude is {apparent_magnitude:.2f}, which is not < {MAGNITUDE_LIMIT}.")
            continue
            
        # If both checks pass
        analysis_log.append(f"{name}: PASSED. DEC={dec} is visible. Apparent magnitude is {apparent_magnitude:.2f} < {MAGNITUDE_LIMIT}.")
        observable_stars.append(name)

    return sorted(observable_stars), analysis_log

def check_final_answer():
    """
    Compares the calculated correct answer with the provided LLM's answer.
    """
    # The final answer from the LLM is 'A', which corresponds to Star3 and Star5.
    llm_answer_option = "A"
    options = {
        "A": ["Star3", "Star5"],
        "B": ["Star4", "Star5"],
        "C": ["Star2", "Star3"],
        "D": ["Star1", "Star4"],
    }
    
    expected_stars = sorted(options.get(llm_answer_option, []))
    
    if not expected_stars:
        return f"Invalid option '{llm_answer_option}' provided."

    calculated_stars, analysis_log = check_star_observability()

    if calculated_stars == expected_stars:
        return "Correct"
    else:
        error_message = f"Incorrect. The provided answer is {expected_stars}, but the code calculated the correct answer to be {calculated_stars}.\n\n"
        error_message += "Here is the step-by-step analysis:\n"
        for log_entry in analysis_log:
            error_message += f"- {log_entry}\n"
        
        # Pinpoint the specific discrepancies
        missing_from_llm = set(calculated_stars) - set(expected_stars)
        extra_in_llm = set(expected_stars) - set(calculated_stars)

        if extra_in_llm:
            error_message += f"\nThe LLM's answer incorrectly includes {list(extra_in_llm)}, which fail the checks.\n"
        if missing_from_llm:
            error_message += f"\nThe LLM's answer fails to include {list(missing_from_llm)}, which pass all checks.\n"
            
        return error_message

# Run the check and print the result
result = check_final_answer()
print(result)