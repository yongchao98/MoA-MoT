import math

def check_correctness():
    """
    Checks the correctness of the LLM's answer by recalculating the velocity of particle A.
    """
    # --- Define Constants ---
    # Rest mass energy of a proton in MeV (using CODATA 2018 value for precision)
    m_p_c2 = 938.27208816
    # Rest mass energy of particle A in MeV, as given in the question
    m_A_c2 = 300.0

    # --- Physics Calculation ---
    # From conservation of energy: 2 * m_p_c2 = 4 * gamma * m_A_c2
    # Solve for the Lorentz factor, gamma
    try:
        gamma = m_p_c2 / (2 * m_A_c2)
    except ZeroDivisionError:
        return "Error: Mass of particle A cannot be zero."

    # The Lorentz factor must be >= 1 for a physical velocity (0 <= v < c)
    if gamma < 1:
        return f"Calculation error: Lorentz factor gamma is {gamma:.4f}, which is less than 1. This implies an impossible physical scenario."

    # Solve for beta (v/c) from the definition of gamma: gamma = 1 / sqrt(1 - beta^2)
    # beta = sqrt(1 - 1/gamma^2)
    beta = math.sqrt(1 - (1 / (gamma**2)))

    # --- Verify the Answer ---
    # The options provided in the question are:
    # A) 0.77c, B) 0.91c, C) 0.96c, D) 0.86c
    options = {
        'A': 0.77,
        'B': 0.91,
        'C': 0.96,
        'D': 0.86
    }
    
    # The final answer given by the LLM is <<<A>>>.
    llm_answer_key = 'A'
    
    # Find which option is numerically closest to our calculated result.
    # This is robust to minor rounding differences.
    closest_option_key = min(options, key=lambda k: abs(options[k] - beta))

    # Check if the LLM's answer matches the correct option based on our calculation.
    if llm_answer_key == closest_option_key:
        return "Correct"
    else:
        return (f"Incorrect. The calculated velocity is v = {beta:.4f}c. "
                f"This is closest to option {closest_option_key} ({options[closest_option_key]}c), "
                f"but the provided answer was {llm_answer_key} ({options[llm_answer_key]}c).")

# Execute the check
result = check_correctness()
print(result)