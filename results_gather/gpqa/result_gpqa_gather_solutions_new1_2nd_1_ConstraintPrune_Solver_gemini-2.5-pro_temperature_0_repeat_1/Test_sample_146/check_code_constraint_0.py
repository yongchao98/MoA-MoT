import math

def check_answer():
    """
    This function checks the correctness of the provided answer to the physics problem.
    It recalculates the velocity of particle A based on the principles of energy conservation
    in special relativity and compares it to the given answer.
    """

    # --- Define Physical Constants and Given Values ---
    # Rest mass energy of a proton in MeV (using a standard value, e.g., from CODATA)
    m_p_c2 = 938.27208816  # in MeV
    # Rest mass energy of particle A in MeV, as given in the question
    m_A_c2 = 300.0  # in MeV

    # --- Step-by-Step Calculation ---

    # 1. Calculate the initial energy (E_initial).
    # The system starts with a proton and an antiproton, assumed to be at rest.
    # Their kinetic energy is negligible ("slowly moving").
    # E_initial = (rest mass energy of proton) + (rest mass energy of antiproton)
    # Since m_p = m_antiproton, E_initial = 2 * m_p_c2
    E_initial = 2 * m_p_c2

    # 2. Formulate the final energy (E_final).
    # The final state has 4 particles (2 A+ and 2 A-), all with the same mass m_A.
    # By conservation of momentum (from an initial state of zero momentum), the energy
    # is distributed equally, so all 4 particles have the same speed 'v' and Lorentz factor 'gamma'.
    # The total energy of a single moving particle is E = gamma * m * c^2.
    # E_final = 4 * (gamma * m_A_c2)
    # E_final = 4 * gamma * 300 MeV = 1200 * gamma MeV

    # 3. Apply conservation of energy: E_initial = E_final
    # 2 * m_p_c2 = 1200 * gamma

    # 4. Solve for the Lorentz factor, gamma.
    try:
        gamma = E_initial / (4 * m_A_c2)
    except ZeroDivisionError:
        return "Calculation Error: The mass of particle A (m_A_c2) cannot be zero."

    # The Lorentz factor must be >= 1.
    if gamma < 1:
        return f"Calculation Error: The calculated Lorentz factor gamma is {gamma:.4f}, which is less than 1 and physically impossible."

    # 5. Solve for the velocity (v) from gamma.
    # The Lorentz factor is defined as gamma = 1 / sqrt(1 - (v/c)^2).
    # Let beta = v/c. Then beta = sqrt(1 - 1/gamma^2).
    try:
        beta_squared = 1 - (1 / (gamma**2))
        calculated_beta = math.sqrt(beta_squared)
    except ValueError:
        return f"Calculation Error: Cannot take the square root of a negative number. This occurred because gamma^2 was less than 1."

    # --- Verify the LLM's Answer ---

    # The final answer provided by the LLM is 'A'.
    llm_answer_letter = 'A'

    # The options as listed in the final LLM response.
    options = {
        'A': 0.77,
        'B': 0.91,
        'C': 0.96,
        'D': 0.86
    }

    # Check if the provided answer letter is a valid option.
    if llm_answer_letter not in options:
        return f"Invalid Answer Format: The provided answer '{llm_answer_letter}' is not a valid option."

    # Get the numerical value corresponding to the LLM's chosen answer.
    llm_answer_value = options[llm_answer_letter]

    # Compare the calculated value with the LLM's answer value.
    # A tolerance is used to account for rounding in the multiple-choice options.
    # A tolerance of 0.005 is appropriate for an answer rounded to two decimal places.
    # For example, our calculated value of ~0.7688 should match 0.77.
    tolerance = 0.005
    if abs(calculated_beta - llm_answer_value) < tolerance:
        return "Correct"
    else:
        # If incorrect, find the option that best matches the calculation.
        best_match_letter = None
        min_diff = float('inf')
        for letter, value in options.items():
            diff = abs(calculated_beta - value)
            if diff < min_diff:
                min_diff = diff
                best_match_letter = letter
        
        reason = (f"Incorrect. The provided answer is '{llm_answer_letter}', which corresponds to a velocity of {llm_answer_value}c. "
                  f"However, the calculation based on energy conservation yields a velocity of approximately {calculated_beta:.4f}c. "
                  f"This calculated value rounds to {round(calculated_beta, 2)}c, which matches option '{best_match_letter}'.")
        return reason

# Run the check and print the result.
print(check_answer())