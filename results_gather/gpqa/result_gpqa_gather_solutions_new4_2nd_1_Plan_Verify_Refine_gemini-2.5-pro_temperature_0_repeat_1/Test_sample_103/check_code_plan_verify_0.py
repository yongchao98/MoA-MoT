import math

def check_exoplanet_period_ratio():
    """
    Checks the correctness of the answer to the exoplanet orbital period question.
    """
    # Given values from the question
    delta_lambda_1 = 5  # miliangstrom
    delta_lambda_2 = 7  # miliangstrom

    # The question asks for the ratio of the orbital period of planet #2 to that of planet #1 (P2 / P1).
    # Based on the physics of the Radial Velocity method for circular orbits and with the given constraints
    # (same star mass, similar planet mass, and assuming similar inclination), the relationship is:
    # P2 / P1 = (Δλ1 / Δλ2)^3
    
    try:
        calculated_ratio = (delta_lambda_1 / delta_lambda_2) ** 3
    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # The provided answer is 'C', which corresponds to the value ~0.36.
    # The options given in the question are:
    # A) ~ 1.40
    # B) ~ 1.96
    # C) ~ 0.36
    # D) ~ 0.85
    
    options = {'A': 1.40, 'B': 1.96, 'C': 0.36, 'D': 0.85}
    llm_choice = 'C'
    
    # Check if the chosen option is the closest to the calculated value.
    closest_option = None
    min_difference = float('inf')

    for option_key, option_value in options.items():
        difference = abs(calculated_ratio - option_value)
        if difference < min_difference:
            min_difference = difference
            closest_option = option_key
            
    if closest_option == llm_choice:
        # The logic in the provided answer is sound, the calculation is correct,
        # and the selected option matches the calculation.
        return "Correct"
    else:
        # The selected option does not match the calculation.
        return (f"Incorrect. The calculated ratio P2/P1 is (5/7)^3 ≈ {calculated_ratio:.4f}. "
                f"The closest option to this value is '{closest_option}' (~{options[closest_option]}), "
                f"but the provided answer was '{llm_choice}' (~{options[llm_choice]}).")

# Execute the check
result = check_exoplanet_period_ratio()
print(result)