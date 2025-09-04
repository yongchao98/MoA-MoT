import math

def check_correctness_of_answer():
    """
    This function checks the correctness of the provided answer to the physics problem.
    It calculates the theoretical values for relative speed and total energy and compares them
    to the values in the selected option.
    """
    # --- Problem Parameters ---
    # Astronaut 1
    m1_m = 2.0  # Mass in units of 'm'
    v1_c = 0.6  # Velocity in units of 'c'
    # Astronaut 2
    m2_m = 3.0  # Mass in units of 'm'
    v2_c = 0.5  # Velocity in units of 'c'

    # --- Theoretical Calculations ---

    # 1. Calculate Relative Speed (v_rel)
    # Using the relativistic velocity subtraction formula: v_rel = (v1 - v2) / (1 - v1*v2/c^2)
    # Since velocities are in units of c, the formula simplifies.
    try:
        calculated_v_rel = (v1_c - v2_c) / (1 - v1_c * v2_c)
    except ZeroDivisionError:
        return "Calculation Error: Division by zero when calculating relative speed."

    # 2. Calculate Total Energy (E_total)
    # The total energy is the sum of individual relativistic energies: E_total = E1 + E2
    # where E = gamma * m_rest * c^2 and gamma = 1 / sqrt(1 - v^2/c^2)
    # We calculate the total energy coefficient for mc^2.
    try:
        # Astronaut 1's energy
        gamma1 = 1 / math.sqrt(1 - v1_c**2)
        energy1_mc2 = gamma1 * m1_m

        # Astronaut 2's energy
        gamma2 = 1 / math.sqrt(1 - v2_c**2)
        energy2_mc2 = gamma2 * m2_m

        # Total system energy
        calculated_e_total = energy1_mc2 + energy2_mc2
    except ValueError:
        return "Calculation Error: Math domain error in Lorentz factor calculation (speed > c?)."
    except ZeroDivisionError:
        return "Calculation Error: Division by zero in Lorentz factor calculation (speed = c?)."

    # --- Compare with the provided answer ---
    # The final answer provided is <<<C>>>, which corresponds to:
    # C) v_rel = 0.14c , E= 5.96 mc^2
    answer_v_rel = 0.14
    answer_e_total = 5.96

    # Check if the calculated values match the answer's values within a tolerance
    # to account for rounding in the options.
    v_rel_is_correct = math.isclose(calculated_v_rel, answer_v_rel, rel_tol=1e-2, abs_tol=1e-3)
    e_total_is_correct = math.isclose(calculated_e_total, answer_e_total, rel_tol=1e-3, abs_tol=1e-3)

    if v_rel_is_correct and e_total_is_correct:
        return "Correct"
    else:
        error_messages = []
        if not v_rel_is_correct:
            error_messages.append(
                f"The relative speed in the answer is incorrect. "
                f"Calculated value: {calculated_v_rel:.4f}c. Answer value: {answer_v_rel}c."
            )
        if not e_total_is_correct:
            error_messages.append(
                f"The total energy in the answer is incorrect. "
                f"Calculated value: {calculated_e_total:.4f}mc^2. Answer value: {answer_e_total}mc^2."
            )
        return "\n".join(error_messages)

# Run the check
print(check_correctness_of_answer())