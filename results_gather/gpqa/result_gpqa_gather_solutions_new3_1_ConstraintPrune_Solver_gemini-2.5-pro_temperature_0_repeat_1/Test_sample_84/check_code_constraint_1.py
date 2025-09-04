import math

def check_answer():
    """
    This function checks the correctness of the provided answer based on the physics principles
    outlined in the problem description.
    """
    # --- Given values from the problem ---
    # Planet 1
    m_p1 = 7.0  # Earth masses
    delta_lambda1 = 0.03  # Angstroms

    # Planet 2
    m_p2 = 5.0  # Earth masses
    delta_lambda2 = 0.04  # Angstroms

    # Options provided in the question
    options = {
        "A": 1.30,
        "B": 1.05,
        "C": 0.98,
        "D": 0.53
    }
    
    # The final answer provided by the LLM to be checked
    llm_answer_letter = "A"
    llm_answer_value = options[llm_answer_letter]

    # --- Calculation ---
    # The ratio of equilibrium temperatures T_eq1 / T_eq2 can be derived as:
    # T_eq1 / T_eq2 = (M_p2 / M_p1) * (K1 / K2)
    # Since the radial velocity semi-amplitude K is proportional to the Doppler shift Δλ,
    # K1 / K2 = Δλ1 / Δλ2.
    # Therefore, T_eq1 / T_eq2 = (M_p2 / M_p1) * (Δλ1 / Δλ2)

    mass_ratio = m_p2 / m_p1
    doppler_shift_ratio = delta_lambda1 / delta_lambda2
    
    calculated_temp_ratio = mass_ratio * doppler_shift_ratio
    
    # --- Verification ---
    # Check if the LLM's answer value is close to the calculated value.
    # We use a tolerance to account for rounding in the options.
    tolerance = 0.01
    
    if abs(llm_answer_value - calculated_temp_ratio) < tolerance:
        return "Correct"
    else:
        # Find the correct option based on the calculation
        correct_option_letter = None
        for letter, value in options.items():
            if abs(value - calculated_temp_ratio) < tolerance:
                correct_option_letter = letter
                break
        
        reason = (
            f"The provided answer is A, which corresponds to a value of {llm_answer_value}.\n"
            f"However, the calculation based on the problem's data yields a different result.\n"
            f"Calculation steps:\n"
            f"1. The ratio of equilibrium temperatures is T_eq1 / T_eq2 = (M_p2 / M_p1) * (Δλ1 / Δλ2).\n"
            f"2. Mass ratio (M_p2 / M_p1) = {m_p2} / {m_p1} = {mass_ratio:.4f}.\n"
            f"3. Doppler shift ratio (Δλ1 / Δλ2) = {delta_lambda1} / {delta_lambda2} = {doppler_shift_ratio:.4f}.\n"
            f"4. Calculated temperature ratio = {mass_ratio:.4f} * {doppler_shift_ratio:.4f} = {calculated_temp_ratio:.4f}.\n"
            f"The calculated value is approximately 0.536, which corresponds to option D (~0.53), not option A (~1.30)."
        )
        return reason

# Run the check
result = check_answer()
print(result)