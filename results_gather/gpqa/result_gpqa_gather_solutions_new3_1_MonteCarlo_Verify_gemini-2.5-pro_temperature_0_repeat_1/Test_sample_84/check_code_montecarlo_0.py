import math

def check_answer():
    """
    This function checks the correctness of the provided answer for the exoplanet temperature ratio problem.
    It recalculates the ratio based on the physical principles and data given in the question.
    """
    # --- Values from the question ---
    # Planet 1
    m_p1 = 7.0  # Earth masses
    delta_lambda1 = 0.03  # Angstroms

    # Planet 2
    m_p2 = 5.0  # Earth masses
    delta_lambda2 = 0.04  # Angstroms

    # --- The final answer provided by the LLM ---
    llm_final_choice = 'B'

    # --- The options given in the question ---
    options = {
        'A': 1.05,
        'B': 0.53,
        'C': 1.30,
        'D': 0.98
    }

    # --- Step-by-step calculation ---
    # The ratio of equilibrium temperatures T_eq1 / T_eq2 can be derived as:
    # T_eq1 / T_eq2 = (M_p2 / M_p1) * (K1 / K2)
    # where K is the radial velocity semi-amplitude.
    # Since K is proportional to the Doppler shift Δλ, we have:
    # T_eq1 / T_eq2 = (M_p2 / M_p1) * (Δλ1 / Δλ2)

    # Calculate the mass ratio
    mass_ratio = m_p2 / m_p1  # This is M_p2 / M_p1

    # Calculate the Doppler shift ratio
    doppler_ratio = delta_lambda1 / delta_lambda2 # This is Δλ1 / Δλ2

    # Calculate the final temperature ratio
    calculated_temp_ratio = mass_ratio * doppler_ratio

    # --- Verification ---
    # Get the numerical value of the LLM's chosen answer
    llm_answer_value = options.get(llm_final_choice)

    if llm_answer_value is None:
        return f"Error: The provided answer choice '{llm_final_choice}' is not one of the valid options (A, B, C, D)."

    # Compare the calculated result with the value of the chosen option.
    # A relative tolerance of 5% is used to account for the '~' in the options.
    if math.isclose(calculated_temp_ratio, llm_answer_value, rel_tol=0.05):
        return "Correct"
    else:
        # Find the correct option based on the calculation
        correct_option = 'None'
        min_diff = float('inf')
        for opt, val in options.items():
            diff = abs(calculated_temp_ratio - val)
            if diff < min_diff:
                min_diff = diff
                correct_option = opt
        
        reason = (
            f"Incorrect. The provided answer is {llm_final_choice} ({llm_answer_value}), "
            f"but the calculated ratio is {calculated_temp_ratio:.4f}.\n"
            f"This value is closest to option {correct_option} ({options[correct_option]}).\n"
            f"Calculation details:\n"
            f"T_eq1 / T_eq2 = (M_p2 / M_p1) * (Δλ1 / Δλ2)\n"
            f"T_eq1 / T_eq2 = ({m_p2} / {m_p1}) * ({delta_lambda1} / {delta_lambda2})\n"
            f"T_eq1 / T_eq2 = {15/28:.4f} ≈ 0.5357"
        )
        return reason

# Run the checker and print the result
result = check_answer()
print(result)