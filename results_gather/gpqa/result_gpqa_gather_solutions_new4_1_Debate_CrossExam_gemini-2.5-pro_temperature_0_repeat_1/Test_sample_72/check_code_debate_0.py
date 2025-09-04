import math

def check_correctness():
    """
    This function checks the correctness of the provided answer to the physics problem.
    It calculates the relative speed and total energy from scratch and compares them
    to the values given in the selected option 'A'.
    """

    # --- Problem Parameters ---
    # Astronaut 1
    m1_coeff = 2.0  # Mass is 2m
    v1_coeff = 0.6  # Velocity is 0.6c

    # Astronaut 2
    m2_coeff = 3.0  # Mass is 3m
    v2_coeff = 0.5  # Velocity is 0.5c

    # --- 1. Calculate the Relative Speed (v_rel) ---
    # The relativistic velocity subtraction formula for objects moving in the same direction is:
    # v_rel = (v1 - v2) / (1 - (v1 * v2) / c^2)
    # Since velocities are given as fractions of c, we can calculate the coefficient of c.
    try:
        calculated_v_rel_coeff = (v1_coeff - v2_coeff) / (1 - v1_coeff * v2_coeff)
    except ZeroDivisionError:
        return "Incorrect. Calculation for relative speed resulted in a division by zero."

    # --- 2. Calculate the Total Energy (E_total) ---
    # The total energy of the system is the sum of the individual relativistic energies.
    # E_total = E1 + E2
    # The formula for a particle's energy is E = gamma * m * c^2, where gamma = 1 / sqrt(1 - (v/c)^2).
    # We will calculate the final coefficient for mc^2.
    try:
        # Energy of Astronaut 1
        gamma1 = 1 / math.sqrt(1 - v1_coeff**2)
        E1_coeff = gamma1 * m1_coeff

        # Energy of Astronaut 2
        gamma2 = 1 / math.sqrt(1 - v2_coeff**2)
        E2_coeff = gamma2 * m2_coeff

        # Total energy coefficient
        calculated_E_total_coeff = E1_coeff + E2_coeff
    except ValueError:
        return "Incorrect. Calculation for gamma factor resulted in a math domain error (e.g., sqrt of a negative number)."


    # --- 3. Extract Values from the Chosen Answer ---
    # The final answer provided is <<<A>>>.
    # The question lists option A as: v_rel = 0.14c, E = 5.96 mc^2
    answer_v_rel_coeff = 0.14
    answer_E_total_coeff = 5.96

    # --- 4. Compare Calculated Values with Answer Values ---
    # A tolerance is used for comparing floating-point numbers, as the options are rounded.
    tolerance = 0.01

    is_v_rel_correct = math.isclose(calculated_v_rel_coeff, answer_v_rel_coeff, rel_tol=tolerance)
    is_E_total_correct = math.isclose(calculated_E_total_coeff, answer_E_total_coeff, rel_tol=tolerance)

    if is_v_rel_correct and is_E_total_correct:
        return "Correct"
    else:
        error_messages = []
        if not is_v_rel_correct:
            error_messages.append(f"The calculated relative speed is approximately {calculated_v_rel_coeff:.3f}c, which does not match the answer's value of {answer_v_rel_coeff}c.")
        if not is_E_total_correct:
            error_messages.append(f"The calculated total energy is approximately {calculated_E_total_coeff:.3f}mc^2, which does not match the answer's value of {answer_E_total_coeff}mc^2.")
        return "Incorrect. " + " ".join(error_messages)

# Execute the check and print the result
result = check_correctness()
print(result)