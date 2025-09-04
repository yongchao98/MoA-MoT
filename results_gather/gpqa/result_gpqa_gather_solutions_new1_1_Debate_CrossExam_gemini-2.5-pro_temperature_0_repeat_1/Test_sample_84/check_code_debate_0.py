import math

def check_planetary_temperature_ratio():
    """
    This function checks the correctness of the provided answer for the planetary temperature ratio problem.
    It calculates the correct ratio based on the given physical parameters and compares it
    to the proposed answer.
    """
    # --- Input Parameters from the Question ---
    # Planet 1
    m_p1 = 7.0  # Earth masses
    delta_lambda1 = 0.03  # Angstroms

    # Planet 2
    m_p2 = 5.0  # Earth masses
    delta_lambda2 = 0.04  # Angstroms

    # The answer to be checked, as provided in the prompt.
    # The prompt's final line states "The final answer is A."
    answer_to_check = 'A'

    # --- Options from the Question ---
    options = {
        "A": 1.05,
        "B": 1.30,
        "C": 0.53,
        "D": 0.98
    }

    # --- Calculation ---
    # The ratio of equilibrium temperatures T_eq1 / T_eq2 can be derived as:
    # T_eq1 / T_eq2 = (M_p2 / M_p1) * (K1 / K2)
    # Since the radial velocity semi-amplitude K is proportional to the Doppler shift Δλ,
    # the formula becomes:
    # T_eq1 / T_eq2 = (M_p2 / M_p1) * (Δλ1 / Δλ2)

    mass_ratio_p2_over_p1 = m_p2 / m_p1
    shift_ratio_1_over_2 = delta_lambda1 / delta_lambda2

    calculated_ratio = mass_ratio_p2_over_p1 * shift_ratio_1_over_2

    # --- Verification ---
    # Find the option that is numerically closest to the calculated result.
    closest_option = None
    min_difference = float('inf')

    for option_key, option_value in options.items():
        difference = abs(calculated_ratio - option_value)
        if difference < min_difference:
            min_difference = difference
            closest_option = option_key

    # Check if the provided answer matches the correct option.
    if answer_to_check == closest_option:
        return "Correct"
    else:
        return (f"Incorrect. The provided answer is {answer_to_check}, which corresponds to a value of {options[answer_to_check]}. "
                f"The calculated ratio is {calculated_ratio:.4f}. "
                f"The closest option is {closest_option} with a value of {options[closest_option]:.2f}.")

# Run the check and print the result.
result = check_planetary_temperature_ratio()
print(result)
