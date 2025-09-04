import math

def check_correctness():
    """
    This function checks the correctness of the selected answer for the given physics problem.
    It calculates the theoretical values for relative speed and total energy and compares them
    against the values from the proposed answer (Option D).
    """
    # --- Problem Parameters ---
    # Astronaut 1: mass = 2m, velocity = 0.6c
    m1_factor = 2.0
    v1_factor = 0.6
    # Astronaut 2: mass = 3m, velocity = 0.5c
    m2_factor = 3.0
    v2_factor = 0.5

    # --- Answer to be Checked (Option D) ---
    # v_rel = 0.14c , E = 5.96 mc^2
    answer_v_rel_factor = 0.14
    answer_e_total_factor = 5.96

    # --- Calculation: Relative Speed ---
    # The relativistic velocity subtraction formula for co-linear motion is:
    # v_rel = (v_a - v_b) / (1 - (v_a * v_b) / c^2)
    # Since velocities are given as factors of c, the formula simplifies.
    # The question asks for speed, so we take the absolute value.
    calculated_v_rel_factor = abs((v1_factor - v2_factor) / (1 - v1_factor * v2_factor))

    # --- Calculation: Total Energy ---
    # The total energy of a system is the sum of the individual relativistic energies.
    # E_total = E1 + E2 = (gamma1 * m1 + gamma2 * m2) * c^2
    # where the Lorentz factor gamma = 1 / sqrt(1 - v^2/c^2)

    # Calculate Lorentz factor for astronaut 1
    gamma1 = 1 / math.sqrt(1 - v1_factor**2)
    # Calculate Lorentz factor for astronaut 2
    gamma2 = 1 / math.sqrt(1 - v2_factor**2)

    # Calculate the total energy factor (the coefficient of mc^2)
    calculated_e_total_factor = (gamma1 * m1_factor) + (gamma2 * m2_factor)

    # --- Verification ---
    # We use math.isclose() to compare the results, allowing for a small tolerance
    # to account for the rounding in the answer options.
    # The calculated v_rel is ~0.1429, which rounds to 0.14. A 3% tolerance is safe.
    # The calculated E_total is ~5.964, which rounds to 5.96. A 1% tolerance is safe.
    v_rel_is_correct = math.isclose(calculated_v_rel_factor, answer_v_rel_factor, rel_tol=0.03)
    e_total_is_correct = math.isclose(calculated_e_total_factor, answer_e_total_factor, rel_tol=0.01)

    if v_rel_is_correct and e_total_is_correct:
        return "Correct"
    else:
        error_messages = []
        if not v_rel_is_correct:
            error_messages.append(
                f"Constraint failed: Relative Speed. "
                f"The answer gives v_rel = {answer_v_rel_factor}c, but the calculated value is approximately {calculated_v_rel_factor:.4f}c."
            )
        if not e_total_is_correct:
            error_messages.append(
                f"Constraint failed: Total Energy. "
                f"The answer gives E_total = {answer_e_total_factor}mc^2, but the calculated value is approximately {calculated_e_total_factor:.4f}mc^2."
            )
        return "\n".join(error_messages)

# Run the check
result = check_correctness()
print(result)