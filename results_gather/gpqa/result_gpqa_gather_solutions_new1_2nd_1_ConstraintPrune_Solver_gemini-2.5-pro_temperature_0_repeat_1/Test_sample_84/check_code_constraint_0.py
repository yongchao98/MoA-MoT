import math

def check_exoplanet_temperature_ratio():
    """
    This function checks the correctness of the LLM's answer by recalculating the ratio of equilibrium temperatures.
    """
    # --- Given values from the question ---
    # Mass of Planet 1 (in Earth masses)
    M_p1 = 7.0
    # Mass of Planet 2 (in Earth masses)
    M_p2 = 5.0
    # Doppler shift for Planet 1 (in Angstroms)
    delta_lambda1 = 0.03
    # Doppler shift for Planet 2 (in Angstroms)
    delta_lambda2 = 0.04

    # --- LLM's final answer and options from the prompt ---
    llm_answer_letter = 'C'
    options = {
        'A': 1.05,
        'B': 1.30,
        'C': 0.53,
        'D': 0.98
    }

    # --- Calculation based on the derived physical formula ---
    # T_eq1 / T_eq2 = (M_p2 / M_p1) * (Δλ1 / Δλ2)
    try:
        calculated_ratio = (M_p2 / M_p1) * (delta_lambda1 / delta_lambda2)
    except ZeroDivisionError:
        return "Error: Division by zero in calculation. Check input values."

    # --- Verification ---
    # Find the option that is numerically closest to the calculated result
    # This handles the "~" (approximately equal) nature of the options
    closest_option_letter = min(options, key=lambda k: abs(options[k] - calculated_ratio))

    # Check if the LLM's chosen letter matches the closest calculated option
    if llm_answer_letter == closest_option_letter:
        # The logic and choice are correct.
        return "Correct"
    else:
        # The LLM chose the wrong option letter.
        llm_choice_value = options.get(llm_answer_letter, "N/A")
        closest_option_value = options[closest_option_letter]
        reason = (f"The answer is incorrect.\n"
                  f"The calculated ratio of equilibrium temperatures (T_eq1 / T_eq2) is {calculated_ratio:.4f}.\n"
                  f"This value is closest to option '{closest_option_letter}' which is ~{closest_option_value}.\n"
                  f"The provided answer was '{llm_answer_letter}' (~{llm_choice_value}), which is not the best match.")
        return reason

# Execute the check
result = check_exoplanet_temperature_ratio()
print(result)