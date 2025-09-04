import math

def check_answer():
    """
    This function checks the correctness of the given answer for the relativistic physics problem.
    """
    # Problem parameters
    v1_c = 0.6  # Velocity of astronaut 1 in units of c
    v2_c = 0.5  # Velocity of astronaut 2 in units of c
    m1_factor = 2  # Mass of astronaut 1 in units of m
    m2_factor = 3  # Mass of astronaut 2 in units of m

    # The options as provided in the question
    options = {
        'A': {'v_rel': 0.14, 'E': 5.96},
        'B': {'v_rel': 0.14, 'E': 5.00},
        'C': {'v_rel': 1.1, 'E': 1.00},
        'D': {'v_rel': 0.1, 'E': 4.96}
    }

    # The final answer provided by the LLM
    llm_answer_choice = 'A'
    
    # --- Step 1: Calculate the true relative speed (v_rel) ---
    # Using the relativistic velocity subtraction formula: v_rel = (v1 - v2) / (1 - v1*v2/c^2)
    try:
        v_rel_numerator = v1_c - v2_c
        v_rel_denominator = 1 - (v1_c * v2_c)
        calculated_v_rel = v_rel_numerator / v_rel_denominator
    except Exception as e:
        return f"An error occurred during relative speed calculation: {e}"

    # --- Step 2: Calculate the true total energy (E_total) ---
    # Using E_total = E1 + E2 = (gamma1 * m1 * c^2) + (gamma2 * m2 * c^2)
    # We calculate the result in units of mc^2
    try:
        # Lorentz factor (gamma): 1 / sqrt(1 - v^2/c^2)
        gamma1 = 1 / math.sqrt(1 - v1_c**2)
        energy1 = gamma1 * m1_factor

        gamma2 = 1 / math.sqrt(1 - v2_c**2)
        energy2 = gamma2 * m2_factor

        calculated_total_energy = energy1 + energy2
    except Exception as e:
        return f"An error occurred during total energy calculation: {e}"

    # --- Step 3: Compare calculated values with the chosen answer ---
    claimed_values = options.get(llm_answer_choice)
    if not claimed_values:
        return f"Invalid answer choice '{llm_answer_choice}'. Please choose from {list(options.keys())}."

    claimed_v_rel = claimed_values['v_rel']
    claimed_E = claimed_values['E']

    # Compare values, rounding to 2 decimal places as in the options
    v_rel_correct = round(calculated_v_rel, 2) == claimed_v_rel
    E_correct = round(calculated_total_energy, 2) == claimed_E

    if v_rel_correct and E_correct:
        return "Correct"
    else:
        error_messages = []
        if not v_rel_correct:
            error_messages.append(
                f"Relative speed is incorrect. Expected: {round(calculated_v_rel, 2)}c, but the answer claims: {claimed_v_rel}c."
            )
        if not E_correct:
            error_messages.append(
                f"Total energy is incorrect. Expected: {round(calculated_total_energy, 2)}mc^2, but the answer claims: {claimed_E}mc^2."
            )
        return "\n".join(error_messages)

# Execute the check
result = check_answer()
print(result)