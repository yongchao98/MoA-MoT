import math

def check_answer_correctness():
    """
    This function checks the correctness of the LLM's answer to the physics problem.
    
    Problem: p + p_bar -> 2A+ + 2A-
    Given: m_A * c^2 = 300 MeV, antiproton is "slowly moving".
    Find: Velocity of particle A.
    Options: A) 0.86c, B) 0.77c, C) 0.91c, D) 0.96c
    LLM's Final Answer to check: <<<B>>>
    """

    # --- Define constants and parameters ---
    
    # Rest energy of a proton in MeV. Using the CODATA 2018 value for accuracy.
    # Minor variations (e.g., 938, 938.3) don't change the result significantly.
    m_p_c2 = 938.27208816  # MeV
    
    # Rest energy of particle A in MeV, as given in the question.
    m_A_c2 = 300.0  # MeV
    
    # The options provided in the question.
    options = {
        'A': 0.86,
        'B': 0.77,
        'C': 0.91,
        'D': 0.96
    }
    
    # The final answer provided by the LLM to be checked.
    llm_final_answer_option = 'B'

    # --- Perform the physics calculation ---

    # 1. Calculate the total initial energy (E_initial).
    # The "slowly moving" antiproton implies negligible initial kinetic energy.
    # The initial energy is the sum of the rest energies of the proton and antiproton.
    E_initial = 2 * m_p_c2

    # 2. Apply conservation of energy (E_initial = E_final) and solve for the Lorentz factor (gamma).
    # The final state has 4 particles of type A. By conservation of momentum
    # from a zero-momentum initial state, all 4 particles must have the same speed 'v'.
    # The total final energy is E_final = 4 * gamma * m_A_c2
    gamma = E_initial / (4 * m_A_c2)
    
    # 3. Calculate the velocity (v) from the Lorentz factor.
    # Let beta = v/c. The relationship is gamma = 1 / sqrt(1 - beta^2).
    # Rearranging gives beta = sqrt(1 - 1/gamma^2).
    if gamma <= 1:
        # This case should not happen if E_initial > 4 * m_A_c2, which it is.
        return f"Incorrect: Calculation resulted in a non-physical Lorentz factor gamma = {gamma:.4f}, which must be greater than 1."
        
    calculated_beta = math.sqrt(1 - (1 / gamma**2))

    # --- Verify the LLM's answer ---

    # Get the numerical value corresponding to the LLM's chosen option.
    if llm_final_answer_option not in options:
        return f"Incorrect: The final answer '{llm_final_answer_option}' is not one of the valid options A, B, C, or D."
        
    llm_answer_value = options[llm_final_answer_option]

    # Compare the calculated result with the LLM's chosen answer.
    # A tolerance is used to account for potential rounding in the option values.
    tolerance = 0.005
    if abs(calculated_beta - llm_answer_value) < tolerance:
        return "Correct"
    else:
        # Find the correct option based on the calculation.
        correct_option = 'None'
        min_diff = float('inf')
        for option, value in options.items():
            diff = abs(calculated_beta - value)
            if diff < min_diff:
                min_diff = diff
                correct_option = option
        
        reason = (f"Incorrect: The calculated velocity is v = {calculated_beta:.4f}c. "
                  f"This value corresponds to option {correct_option} ({options.get(correct_option)}c). "
                  f"The provided answer was option {llm_final_answer_option} ({llm_answer_value}c).")
        
        return reason

# Execute the check and print the result.
result = check_answer_correctness()
print(result)