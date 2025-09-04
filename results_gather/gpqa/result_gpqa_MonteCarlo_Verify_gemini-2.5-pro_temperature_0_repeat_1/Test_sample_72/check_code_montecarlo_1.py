import math

def check_relativity_answer():
    """
    This function calculates the relative speed and total energy for the given physics problem
    and compares the results with the values from the proposed answer (Option A).
    """
    # --- Define problem parameters from the question ---
    # Astronaut 1: mass 2m, speed 0.6c
    m1_factor = 2.0
    v1 = 0.6  # Speed in units of c

    # Astronaut 2: mass 3m, speed 0.5c
    m2_factor = 3.0
    v2 = 0.5  # Speed in units of c

    # --- Values from the proposed answer (Option A) ---
    v_rel_answer = 0.14
    E_total_answer = 5.96

    # --- 1. Calculate the theoretical relative speed ---
    # The relativistic velocity subtraction formula for collinear motion is:
    # v_rel = (v_a - v_b) / (1 - (v_a * v_b) / c^2)
    # Since speeds are in units of c, we can set c=1.
    try:
        v_rel_calculated = (v1 - v2) / (1 - v1 * v2)
    except ZeroDivisionError:
        return "Calculation Error: Division by zero in relative velocity formula."

    # --- 2. Calculate the theoretical total energy ---
    # The total energy of a particle is E = γ * m_rest * c^2, where γ = 1 / sqrt(1 - v^2/c^2).
    # The total energy of the system is the sum of the individual energies.
    # We calculate the energy in units of mc^2.
    try:
        # Energy of Astronaut 1
        gamma1 = 1 / math.sqrt(1 - v1**2)
        E1 = gamma1 * m1_factor

        # Energy of Astronaut 2
        gamma2 = 1 / math.sqrt(1 - v2**2)
        E2 = gamma2 * m2_factor

        # Total system energy
        E_total_calculated = E1 + E2
    except ValueError:
        return "Calculation Error: Cannot take the square root of a negative number (speed >= c)."

    # --- 3. Compare calculated values with the answer ---
    # A tolerance is necessary because the values in the options are rounded.
    # A tolerance of 0.01 is appropriate for the given precision.
    tolerance = 0.01

    is_v_rel_correct = math.isclose(v_rel_calculated, v_rel_answer, abs_tol=tolerance)
    is_E_total_correct = math.isclose(E_total_calculated, E_total_answer, abs_tol=tolerance)

    # --- 4. Return the result of the check ---
    if is_v_rel_correct and is_E_total_correct:
        return "Correct"
    else:
        error_messages = []
        if not is_v_rel_correct:
            error_messages.append(
                f"The relative speed is incorrect. The answer states v_rel = {v_rel_answer}c, but the calculated value is approximately {v_rel_calculated:.4f}c."
            )
        if not is_E_total_correct:
            error_messages.append(
                f"The total energy is incorrect. The answer states E = {E_total_answer} mc^2, but the calculated value is approximately {E_total_calculated:.4f} mc^2."
            )
        return "\n".join(error_messages)

# Execute the check and print the result
result = check_relativity_answer()
print(result)