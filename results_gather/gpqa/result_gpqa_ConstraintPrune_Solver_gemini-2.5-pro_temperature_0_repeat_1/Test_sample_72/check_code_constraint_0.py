import math

def check_correctness():
    """
    This function checks the correctness of the given answer to the physics problem.
    It calculates the theoretical values for relative speed and total energy and
    compares them against the values provided in the selected option D.
    """
    # --- Given parameters from the question ---
    # We can treat m and c as units, so we work with the numerical factors.
    m1_factor = 2.0
    v1_factor = 0.6  # in units of c
    m2_factor = 3.0
    v2_factor = 0.5  # in units of c

    # --- Values from the proposed answer (Option D) ---
    answer_v_rel_factor = 0.14
    answer_E_total_factor = 5.96

    # --- Step 1: Calculate the theoretical relative speed ---
    # The relativistic velocity subtraction formula for collinear motion is:
    # v_rel = (v1 - v2) / (1 - (v1 * v2) / c^2)
    # In terms of factors of c:
    try:
        calculated_v_rel_factor = (v1_factor - v2_factor) / (1 - v1_factor * v2_factor)
    except ZeroDivisionError:
        return "Calculation Error: Division by zero when calculating relative velocity."

    # --- Step 2: Calculate the theoretical total energy ---
    # The total energy of a system is the sum of the individual relativistic energies.
    # E_total = E1 + E2 = (gamma1 * m1 + gamma2 * m2) * c^2
    # where gamma = 1 / sqrt(1 - v^2/c^2)
    try:
        gamma1 = 1 / math.sqrt(1 - v1_factor**2)
        gamma2 = 1 / math.sqrt(1 - v2_factor**2)
        calculated_E_total_factor = (gamma1 * m1_factor) + (gamma2 * m2_factor)
    except ValueError:
        return "Calculation Error: A velocity is greater than or equal to c, cannot compute Lorentz factor."

    # --- Step 3: Compare calculated values with the answer ---
    # We use math.isclose() to account for potential rounding in the options.
    # The options are given to 2 decimal places, so a relative tolerance of 1% is reasonable.
    v_rel_is_correct = math.isclose(answer_v_rel_factor, calculated_v_rel_factor, rel_tol=0.01)
    E_total_is_correct = math.isclose(answer_E_total_factor, calculated_E_total_factor, rel_tol=0.01)

    if v_rel_is_correct and E_total_is_correct:
        return "Correct"
    else:
        error_messages = []
        if not v_rel_is_correct:
            error_messages.append(
                f"Constraint failed: Relative speed.\n"
                f"The answer states v_rel = {answer_v_rel_factor}c, but the calculated value is approximately {calculated_v_rel_factor:.4f}c."
            )
        if not E_total_is_correct:
            error_messages.append(
                f"Constraint failed: Total energy.\n"
                f"The answer states E = {answer_E_total_factor}mc^2, but the calculated value is approximately {calculated_E_total_factor:.4f}mc^2."
            )
        return "\n".join(error_messages)

# Run the check and print the result
result = check_correctness()
print(result)