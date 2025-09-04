import math

def check_physics_problem():
    """
    This function checks the correctness of the answer to the given relativistic physics problem.

    The problem:
    Two astronauts of masses 2m and 3m travel at relativistic speeds of 0.6c and 0.5c,
    respectively, in the positive x-direction. What will be the relative speed and
    the total energy of the system?

    The proposed answer (Option C):
    v_rel = 0.14c, E = 5.96 mc^2
    """

    # --- Define constants and given values ---
    # We can work in units where m=1 and c=1, as the final answer is expressed in terms of m and c.
    m = 1.0
    c = 1.0
    
    m1 = 2 * m
    v1 = 0.6 * c
    
    m2 = 3 * m
    v2 = 0.5 * c

    # --- 1. Calculate the relative speed (v_rel) ---
    # The relativistic velocity subtraction formula for collinear motion is:
    # v_rel = (u - v) / (1 - uv/c^2)
    # Here, u = v1 and v = v2.
    try:
        v_rel_numerator = v1 - v2
        v_rel_denominator = 1 - (v1 * v2) / (c**2)
        calculated_v_rel = v_rel_numerator / v_rel_denominator
    except Exception as e:
        return f"An error occurred during relative speed calculation: {e}"

    # --- 2. Calculate the total energy of the system (E_total) ---
    # The total energy of a system is the sum of the individual relativistic energies.
    # E_particle = gamma * mass_rest * c^2
    # where gamma = 1 / sqrt(1 - v^2/c^2)
    try:
        # Energy of the first astronaut
        gamma1 = 1 / math.sqrt(1 - (v1/c)**2)
        E1 = gamma1 * m1 * c**2

        # Energy of the second astronaut
        gamma2 = 1 / math.sqrt(1 - (v2/c)**2)
        E2 = gamma2 * m2 * c**2

        # Total energy of the system (in units of mc^2)
        calculated_E_total = (E1 + E2) / (m * c**2)
    except Exception as e:
        return f"An error occurred during total energy calculation: {e}"

    # --- 3. Compare calculated values with the given answer (Option C) ---
    answer_v_rel = 0.14   # in units of c
    answer_E_total = 5.96 # in units of mc^2

    # Use a tolerance for floating-point comparison, as the answer is rounded.
    # A tolerance of 1% (rel_tol=0.01) is reasonable for values rounded to two decimal places.
    is_v_rel_correct = math.isclose(calculated_v_rel, answer_v_rel, rel_tol=0.01)
    is_E_total_correct = math.isclose(calculated_E_total, answer_E_total, rel_tol=0.01)

    # --- 4. Formulate the final result ---
    if is_v_rel_correct and is_E_total_correct:
        return "Correct"
    else:
        error_messages = []
        if not is_v_rel_correct:
            error_messages.append(
                f"Relative speed is incorrect. "
                f"Calculated v_rel = {calculated_v_rel:.4f}c, but the answer states {answer_v_rel}c."
            )
        if not is_E_total_correct:
            error_messages.append(
                f"Total energy is incorrect. "
                f"Calculated E_total = {calculated_E_total:.4f}mc^2, but the answer states {answer_E_total}mc^2."
            )
        return "\n".join(error_messages)

# Execute the check
result = check_physics_problem()
print(result)