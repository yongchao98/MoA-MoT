import math

def check_answer():
    """
    This function checks the correctness of the given answer to the relativistic physics problem.
    It calculates the values based on the problem statement and compares them to the provided answer.
    """
    # --- Problem Parameters ---
    # Masses are given as multiples of a base mass 'm'
    m1_factor = 2.0
    m2_factor = 3.0
    # Velocities are given as fractions of the speed of light 'c'
    v1_factor = 0.6
    v2_factor = 0.5

    # --- The Answer to be Checked (Option D) ---
    # v_rel = 0.14c , E= 5.96 mc^2
    expected_v_rel_factor = 0.14
    expected_E_total_factor = 5.96

    # --- Independent Calculation ---

    # 1. Calculate the relative speed (v_rel)
    # The relativistic velocity subtraction formula is v_rel = (v2 - v1) / (1 - v1*v2/c^2).
    # Since velocities are given as factors of c, the formula simplifies.
    # The question asks for speed, which is the magnitude of the velocity.
    try:
        calculated_v_rel_factor = abs((v2_factor - v1_factor) / (1 - v1_factor * v2_factor))
    except ZeroDivisionError:
        return "Calculation failed: Division by zero in relative speed formula. This occurs if v1*v2 = c^2."


    # 2. Calculate the total energy (E_total)
    # The relativistic energy formula is E = gamma * m0 * c^2, where gamma = 1 / sqrt(1 - v^2/c^2).
    # The total energy of the system is the sum of the individual relativistic energies.
    # E_total = (gamma1 * m1 + gamma2 * m2) * c^2
    try:
        # Calculate Lorentz factor (gamma) for astronaut 1
        gamma1 = 1 / math.sqrt(1 - v1_factor**2)
        # Calculate Lorentz factor (gamma) for astronaut 2
        gamma2 = 1 / math.sqrt(1 - v2_factor**2)
    except ValueError:
        return "Calculation failed: Cannot take square root of a negative number. This occurs if a speed is > c."

    # The total energy factor is the multiplier for mc^2
    calculated_E_total_factor = (gamma1 * m1_factor) + (gamma2 * m2_factor)

    # --- Verification ---
    # We check if the calculated values, when rounded to two decimal places, match the answer.
    # This is a robust way to handle potential floating-point inaccuracies and rounding in the provided options.
    v_rel_matches = round(calculated_v_rel_factor, 2) == expected_v_rel_factor
    E_total_matches = round(calculated_E_total_factor, 2) == expected_E_total_factor

    if v_rel_matches and E_total_matches:
        return "Correct"
    else:
        error_messages = []
        if not v_rel_matches:
            error_messages.append(
                f"Constraint for relative speed is not satisfied. "
                f"Calculated value: {calculated_v_rel_factor:.4f}c, which rounds to {round(calculated_v_rel_factor, 2)}c. "
                f"Answer provided: {expected_v_rel_factor}c."
            )
        if not E_total_matches:
            error_messages.append(
                f"Constraint for total energy is not satisfied. "
                f"Calculated value: {calculated_E_total_factor:.4f}mc^2, which rounds to {round(calculated_E_total_factor, 2)}mc^2. "
                f"Answer provided: {expected_E_total_factor}mc^2."
            )
        return "\n".join(error_messages)

# Run the check
result = check_answer()
print(result)