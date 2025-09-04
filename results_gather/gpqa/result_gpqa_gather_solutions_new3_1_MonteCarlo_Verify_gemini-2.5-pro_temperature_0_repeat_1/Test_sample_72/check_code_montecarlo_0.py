import math

def check_answer():
    """
    This function checks the correctness of the provided answer to the physics problem.
    """
    # Given parameters from the question
    m1_factor = 2.0
    v1_c = 0.6  # v1 in terms of c
    m2_factor = 3.0
    v2_c = 0.5  # v2 in terms of c

    # The final answer provided by the LLM is B
    # Values from option B: v_rel=0.14c, E=5.96mc^2
    proposed_v_rel_c = 0.14
    proposed_E_total_mc2 = 5.96

    # --- Step 1: Calculate the correct relative speed ---
    # The relativistic velocity subtraction formula for co-linear motion is:
    # v_rel = (v1 - v2) / (1 - (v1 * v2) / c^2)
    # In terms of fractions of c, this simplifies to:
    # v_rel_c = (v1_c - v2_c) / (1 - v1_c * v2_c)
    
    # Ensure the formula uses the difference in velocities, the order doesn't matter for the magnitude.
    # We can use abs() to be safe, but v1 > v2 so it's not strictly necessary here.
    calculated_v_rel_c = (v1_c - v2_c) / (1 - v1_c * v2_c)

    # --- Step 2: Calculate the correct total energy ---
    # The total energy of a relativistic particle is E = gamma * m_rest * c^2
    # where gamma = 1 / sqrt(1 - (v/c)^2)
    # The total energy of the system is the sum of the individual energies:
    # E_total = E1 + E2 = (gamma1 * m1 * c^2) + (gamma2 * m2 * c^2)
    # We calculate the total energy in units of mc^2:
    # E_total / (m * c^2) = (gamma1 * m1_factor) + (gamma2 * m2_factor)

    gamma1 = 1 / math.sqrt(1 - v1_c**2)
    gamma2 = 1 / math.sqrt(1 - v2_c**2)

    calculated_E_total_mc2 = (gamma1 * m1_factor) + (gamma2 * m2_factor)

    # --- Step 3: Compare the calculated values with the proposed answer ---
    # We use a tolerance because the options are rounded to two decimal places.
    # A tolerance of 0.01 is appropriate.
    tolerance = 0.01

    v_rel_is_correct = math.isclose(calculated_v_rel_c, proposed_v_rel_c, abs_tol=tolerance)
    E_total_is_correct = math.isclose(calculated_E_total_mc2, proposed_E_total_mc2, abs_tol=tolerance)

    if v_rel_is_correct and E_total_is_correct:
        return "Correct"
    else:
        error_messages = []
        if not v_rel_is_correct:
            error_messages.append(
                f"The relative speed is incorrect. "
                f"The calculated value is v_rel ≈ {calculated_v_rel_c:.4f}c, "
                f"which rounds to {round(calculated_v_rel_c, 2)}c. "
                f"The answer states v_rel = {proposed_v_rel_c}c."
            )
        if not E_total_is_correct:
            error_messages.append(
                f"The total energy is incorrect. "
                f"The calculated value is E ≈ {calculated_E_total_mc2:.4f}mc^2, "
                f"which rounds to {round(calculated_E_total_mc2, 2)}mc^2. "
                f"The answer states E = {proposed_E_total_mc2}mc^2."
            )
        return "\n".join(error_messages)

# Run the check and print the result
result = check_answer()
print(result)