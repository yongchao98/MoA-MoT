import math

def check_correctness_of_answer():
    """
    Checks the correctness of the LLM's answer to the physics problem by
    re-calculating the result from first principles.

    The problem is:
    p + p_bar -> 2A+ + 2A-
    m_A*c^2 = 300 MeV
    Find the velocity of particle A.

    The LLM's final answer is <<<C>>>, which corresponds to 0.77c.
    """

    # --- Problem Constants and Given Values ---
    # Rest mass energy of a proton in MeV (CODATA 2018 value for higher precision)
    m_p_c2 = 938.27208816
    # Rest mass energy of particle A in MeV (from the question)
    m_A_c2 = 300.0

    # --- Multiple Choice Options from the Question ---
    # The question text has: A) 0.86c, B) 0.91c, C) 0.77c, D) 0.96c
    options = {
        'A': 0.86,
        'B': 0.91,
        'C': 0.77,
        'D': 0.96
    }

    # --- LLM's Answer to Verify ---
    # The final answer provided in the prompt is <<<C>>>
    llm_choice = 'C'

    # --- Physics Calculation from First Principles ---
    # 1. Apply conservation of energy: E_initial = E_final
    # E_initial = rest energy of proton + rest energy of antiproton (assuming negligible initial KE)
    e_initial = 2 * m_p_c2

    # E_final = 4 * (total energy of one A particle) = 4 * gamma * m_A_c2
    # So, 2 * m_p_c2 = 4 * gamma * m_A_c2

    # 2. Solve for the Lorentz factor, gamma
    # gamma = (2 * m_p_c2) / (4 * m_A_c2) = m_p_c2 / (2 * m_A_c2)
    try:
        gamma = m_p_c2 / (2 * m_A_c2)
    except ZeroDivisionError:
        return "Error in calculation: Division by zero. m_A_c2 cannot be zero."

    # 3. Check if the process is physically possible (gamma must be >= 1)
    if gamma < 1:
        return (f"Calculation error: Lorentz factor gamma ({gamma:.4f}) is less than 1. "
                f"This implies the total rest mass of the products ({4 * m_A_c2} MeV) "
                f"is greater than the total initial energy ({e_initial:.4f} MeV), "
                f"which is physically impossible.")

    # 4. Solve for the velocity (beta = v/c) from gamma
    # gamma = 1 / sqrt(1 - beta^2)  =>  beta = sqrt(1 - 1/gamma^2)
    try:
        beta_squared = 1 - (1 / (gamma**2))
        calculated_beta = math.sqrt(beta_squared)
    except ValueError:
        return "Calculation error: Cannot take the square root of a negative number while calculating beta."

    # --- Verification Step ---
    # Find which of the given options is numerically closest to the calculated result.
    try:
        closest_option_key = min(options.keys(), key=lambda k: abs(options[k] - calculated_beta))
    except (ValueError, TypeError):
        return "Error in verification: Could not compare calculated value with options."

    # Compare the LLM's choice with the calculated best choice.
    if llm_choice == closest_option_key:
        return "Correct"
    else:
        reason = (
            f"The provided answer choice '{llm_choice}' is incorrect.\n"
            f"The calculated velocity is v/c = {calculated_beta:.4f}.\n"
            f"The value corresponding to the chosen answer '{llm_choice}' is {options[llm_choice]}c.\n"
            f"The option closest to the calculated value is '{closest_option_key}', which corresponds to a velocity of {options[closest_option_key]}c.\n"
            f"\nCalculation details:\n"
            f"1. Initial energy E_i = 2 * m_p*c^2 = {e_initial:.3f} MeV.\n"
            f"2. Final energy E_f = 4 * gamma * m_A*c^2 = 1200 * gamma MeV.\n"
            f"3. Equating energies gives gamma = E_i / 1200 = {gamma:.4f}.\n"
            f"4. Solving for velocity gives v/c = sqrt(1 - 1/gamma^2) = {calculated_beta:.4f}."
        )
        return reason

# print(check_correctness_of_answer())