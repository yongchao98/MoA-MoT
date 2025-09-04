import math

def check_correctness_of_astro_answer():
    """
    This function checks the correctness of the final answer provided for the exoplanet temperature ratio question.
    It recalculates the ratio from the given parameters and compares it to the value and option selected in the final answer.
    """
    
    # --- Given parameters from the question ---
    # Mass of Planet 1 in Earth masses
    M_p1 = 7
    # Mass of Planet 2 in Earth masses
    M_p2 = 5
    # Doppler shift for Planet 1 in Angstroms
    delta_lambda1 = 0.03
    # Doppler shift for Planet 2 in Angstroms
    delta_lambda2 = 0.04

    # --- Information from the final answer to be checked ---
    # The final answer concludes the result is ~0.536 and selects option D.
    final_answer_option = 'D'
    # The options provided in the question prompt
    options = {
        'A': 1.30,
        'B': 1.05,
        'C': 0.98,
        'D': 0.53
    }

    # --- Calculation based on physics principles ---
    # The ratio of equilibrium temperatures T_eq1 / T_eq2 can be derived as:
    # (M_p2 / M_p1) * (delta_lambda1 / delta_lambda2)
    
    mass_ratio = M_p2 / M_p1
    doppler_shift_ratio = delta_lambda1 / delta_lambda2
    
    calculated_ratio = mass_ratio * doppler_shift_ratio
    
    # --- Verification ---
    # 1. Check if the selected option is valid
    if final_answer_option not in options:
        return f"Invalid option: The final answer selected '{final_answer_option}', which is not one of the available options (A, B, C, D)."

    # 2. Check if the calculated numerical value matches the value of the selected option
    expected_value_for_option = options[final_answer_option]
    
    # We use a tolerance because the options are given with '~' (approximately)
    # A tolerance of 0.01 is reasonable for values given to two decimal places.
    if math.isclose(calculated_ratio, expected_value_for_option, abs_tol=0.01):
        return "Correct"
    else:
        # Find which option the calculation actually matches
        correct_option = 'None'
        for opt, val in options.items():
            if math.isclose(calculated_ratio, val, abs_tol=0.01):
                correct_option = opt
                break
        
        return (f"Incorrect. The calculation is correct, yielding a ratio of {calculated_ratio:.4f}. "
                f"This value corresponds to option {correct_option} (~{options.get(correct_option, 'N/A')}). "
                f"However, the final answer selected option {final_answer_option} (~{expected_value_for_option}).")

# Execute the check
result = check_correctness_of_astro_answer()
print(result)