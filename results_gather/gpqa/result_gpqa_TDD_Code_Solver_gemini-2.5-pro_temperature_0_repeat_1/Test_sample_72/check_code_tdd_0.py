import math

def check_answer():
    """
    This function checks the correctness of the given answer to the relativistic physics problem.
    It calculates the relative speed and total energy from first principles and compares them
    to the values in option D.
    """
    # --- Problem Parameters ---
    # Masses are given as multiples of a base mass 'm'
    m1_factor = 2.0
    m2_factor = 3.0
    # Velocities are given as fractions of the speed of light 'c'
    v1_factor = 0.6
    v2_factor = 0.5

    # --- Answer to be Checked (Option D) ---
    expected_v_rel_factor = 0.14
    expected_E_total_factor = 5.96

    # --- Calculation ---

    # 1. Calculate the relative speed (v_rel)
    # The relativistic velocity subtraction formula is: v_rel = (v2 - v1) / (1 - (v1*v2)/c^2)
    # Since velocities are factors of c, the formula simplifies to:
    # v_rel_factor = (v2_factor - v1_factor) / (1 - v1_factor * v2_factor)
    # The question asks for speed, which is the magnitude (absolute value).
    try:
        calculated_v_rel_factor = abs((v2_factor - v1_factor) / (1 - v1_factor * v2_factor))
    except ZeroDivisionError:
        return "Calculation failed: Division by zero in relative speed formula."

    # 2. Calculate the total energy (E_total)
    # The total energy of a particle is E = gamma * m0 * c^2, where gamma = 1 / sqrt(1 - v^2/c^2).
    # The total energy of the system is the sum of the individual energies:
    # E_total = E1 + E2 = (gamma1 * m1 + gamma2 * m2) * c^2
    # We will calculate the factor that multiplies 'mc^2'.

    # Calculate Lorentz factor (gamma) for astronaut 1
    try:
        gamma1 = 1 / math.sqrt(1 - v1_factor**2)
    except ValueError:
        return f"Calculation failed: Invalid velocity v1={v1_factor}c (must be < c)."

    # Calculate Lorentz factor (gamma) for astronaut 2
    try:
        gamma2 = 1 / math.sqrt(1 - v2_factor**2)
    except ValueError:
        return f"Calculation failed: Invalid velocity v2={v2_factor}c (must be < c)."

    # Calculate the total energy factor
    calculated_E_total_factor = (gamma1 * m1_factor) + (gamma2 * m2_factor)

    # --- Verification ---
    # We use math.isclose() to compare floating-point numbers, allowing for small rounding differences.
    # The tolerance is set based on the precision of the given answer (2 decimal places).
    v_rel_is_correct = math.isclose(calculated_v_rel_factor, expected_v_rel_factor, rel_tol=1e-2, abs_tol=0.005)
    E_total_is_correct = math.isclose(calculated_E_total_factor, expected_E_total_factor, rel_tol=1e-3, abs_tol=0.005)

    if v_rel_is_correct and E_total_is_correct:
        return "Correct"
    else:
        error_messages = []
        if not v_rel_is_correct:
            error_messages.append(
                f"Relative speed is incorrect. The answer states v_rel = {expected_v_rel_factor}c, "
                f"but the calculated value is approximately {calculated_v_rel_factor:.4f}c."
            )
        if not E_total_is_correct:
            error_messages.append(
                f"Total energy is incorrect. The answer states E = {expected_E_total_factor}mc^2, "
                f"but the calculated value is approximately {calculated_E_total_factor:.4f}mc^2."
            )
        return "\n".join(error_messages)

# Run the check and print the result
result = check_answer()
print(result)