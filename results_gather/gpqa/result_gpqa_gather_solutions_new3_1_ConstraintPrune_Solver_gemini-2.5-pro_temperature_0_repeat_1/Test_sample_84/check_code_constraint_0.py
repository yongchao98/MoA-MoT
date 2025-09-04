import math

def check_planetary_temperature_ratio():
    """
    This function verifies the solution to the exoplanet temperature ratio problem.
    It calculates the correct ratio based on the provided physical parameters and checks
    which of the multiple-choice options is correct.
    """
    # --- Given Data ---
    # Planet 1
    m_p1 = 7.0  # Earth masses
    delta_lambda1 = 0.03  # Angstroms

    # Planet 2
    m_p2 = 5.0  # Earth masses
    delta_lambda2 = 0.04  # Angstroms

    # --- Physics Derivation ---
    # The ratio of equilibrium temperatures T_eq1 / T_eq2 is sqrt(a2 / a1).
    # The radial velocity semi-amplitude K is proportional to M_p / sqrt(a).
    # This implies sqrt(a) is proportional to M_p / K.
    # Therefore, T_eq1 / T_eq2 = (M_p2 / K2) / (M_p1 / K1) = (M_p2 / M_p1) * (K1 / K2).
    # Since K is proportional to the Doppler shift delta_lambda, K1/K2 = delta_lambda1/delta_lambda2.
    # So, the final formula is:
    # T_eq1 / T_eq2 = (m_p2 / m_p1) * (delta_lambda1 / delta_lambda2)

    # --- Calculation ---
    mass_ratio = m_p2 / m_p1
    doppler_shift_ratio = delta_lambda1 / delta_lambda2
    temperature_ratio = mass_ratio * doppler_shift_ratio

    # --- Compare with Options ---
    options = {
        'A': 1.30,
        'B': 1.05,
        'C': 0.98,
        'D': 0.53
    }

    # Find the closest option to the calculated result
    closest_option = None
    min_difference = float('inf')
    for option_key, option_value in options.items():
        difference = abs(temperature_ratio - option_value)
        if difference < min_difference:
            min_difference = difference
            closest_option = option_key
            
    # The provided answer from the LLMs is a mix, but let's check against the most common incorrect answer 'A'
    # and the correct answer 'D'.
    llm_answer_to_check = 'A' # This is one of the incorrect answers provided.

    if llm_answer_to_check == closest_option:
        return "Correct"
    else:
        error_message = (
            f"Incorrect. The provided answer is {llm_answer_to_check}, but the calculated correct answer is {closest_option}.\n"
            f"Reasoning:\n"
            f"1. The ratio of equilibrium temperatures T_eq1/T_eq2 simplifies to sqrt(a2/a1), where 'a' is the semi-major axis.\n"
            f"2. The radial velocity semi-amplitude K is proportional to M_p/sqrt(a). This leads to the relationship T_eq1/T_eq2 = (M_p2/M_p1) * (K1/K2).\n"
            f"3. The ratio of radial velocities K1/K2 is equal to the ratio of Doppler shifts Δλ1/Δλ2.\n"
            f"4. Calculation: T_eq1/T_eq2 = (5/7) * (0.03/0.04) = (5/7) * (3/4) = 15/28 ≈ {temperature_ratio:.4f}.\n"
            f"5. The calculated value ~0.536 is closest to option D (0.53), not option {llm_answer_to_check} ({options[llm_answer_to_check]})."
        )
        return error_message

# Execute the check and print the result
result = check_planetary_temperature_ratio()
print(result)
