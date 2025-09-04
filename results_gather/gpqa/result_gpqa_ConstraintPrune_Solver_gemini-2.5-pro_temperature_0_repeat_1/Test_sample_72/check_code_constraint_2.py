import math

def check_physics_answer():
    """
    This function checks the correctness of the given answer to the relativistic physics problem.
    """
    # --- Given values from the question ---
    # We work with factors of m and c.
    m1_factor = 2.0  # Mass of astronaut 1 is 2m
    v1_factor = 0.6  # Velocity of astronaut 1 is 0.6c
    m2_factor = 3.0  # Mass of astronaut 2 is 3m
    v2_factor = 0.5  # Velocity of astronaut 2 is 0.5c

    # --- Answer to be checked (from Option D) ---
    answer_v_rel_factor = 0.14
    answer_E_total_factor = 5.96

    # --- Constraint 1: Calculate Relative Velocity ---
    # The formula for relativistic velocity subtraction for objects moving in the same direction is:
    # v_rel = (v1 - v2) / (1 - (v1 * v2) / c^2)
    # We can calculate the factor by which c is multiplied.
    try:
        calc_v_rel_factor = (v1_factor - v2_factor) / (1 - v1_factor * v2_factor)
    except ZeroDivisionError:
        return "Error: Division by zero in relative velocity calculation."

    # --- Constraint 2: Calculate Total Energy ---
    # The total energy of the system is the sum of the individual relativistic energies.
    # E_total = E1 + E2
    # E = gamma * m * c^2, where gamma = 1 / sqrt(1 - v^2/c^2)
    
    # Calculate gamma and energy for astronaut 1
    try:
        gamma1 = 1 / math.sqrt(1 - v1_factor**2)
        E1_factor = gamma1 * m1_factor  # Energy in units of mc^2
    except ValueError:
        return "Error: Math domain error for astronaut 1 (velocity likely >= c)."

    # Calculate gamma and energy for astronaut 2
    try:
        gamma2 = 1 / math.sqrt(1 - v2_factor**2)
        E2_factor = gamma2 * m2_factor  # Energy in units of mc^2
    except ValueError:
        return "Error: Math domain error for astronaut 2 (velocity likely >= c)."

    # Calculate total energy factor
    calc_E_total_factor = E1_factor + E2_factor

    # --- Verification ---
    # The options are rounded to two decimal places, so we round our calculated values
    # to the same precision for a fair comparison.
    rounded_calc_v_rel = round(calc_v_rel_factor, 2)
    rounded_calc_E_total = round(calc_E_total_factor, 2)

    # Check if the rounded calculated values match the answer's values
    v_rel_is_correct = (rounded_calc_v_rel == answer_v_rel_factor)
    E_total_is_correct = (rounded_calc_E_total == answer_E_total_factor)

    if v_rel_is_correct and E_total_is_correct:
        return "Correct"
    else:
        error_messages = []
        if not v_rel_is_correct:
            error_messages.append(
                f"Relative velocity constraint is not satisfied. "
                f"Calculated value: {calc_v_rel_factor:.4f}c (rounds to {rounded_calc_v_rel}c), "
                f"Answer value: {answer_v_rel_factor}c."
            )
        if not E_total_is_correct:
            error_messages.append(
                f"Total energy constraint is not satisfied. "
                f"Calculated value: {calc_E_total_factor:.4f}mc^2 (rounds to {rounded_calc_E_total}mc^2), "
                f"Answer value: {answer_E_total_factor}mc^2."
            )
        return "\n".join(error_messages)

# Run the check
result = check_physics_answer()
print(result)