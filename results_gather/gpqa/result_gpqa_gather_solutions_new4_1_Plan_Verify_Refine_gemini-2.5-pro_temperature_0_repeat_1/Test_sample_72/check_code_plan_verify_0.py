import math

def check_correctness_of_answer():
    """
    This function checks the correctness of the given answer to the physics problem.
    It recalculates the relative speed and total energy from scratch and compares them
    to the values in the selected option.
    """

    # --- Problem Parameters ---
    # Astronaut 1
    m1_coeff = 2.0
    v1_c = 0.6  # Velocity as a fraction of c

    # Astronaut 2
    m2_coeff = 3.0
    v2_c = 0.5  # Velocity as a fraction of c

    # --- Step 1: Calculate the theoretical relative speed ---
    # The relativistic velocity subtraction formula for objects moving in the same direction is:
    # v_rel = (v1 - v2) / (1 - (v1 * v2) / c^2)
    # We calculate the coefficient of c.
    try:
        v_rel_calc_coeff = (v1_c - v2_c) / (1 - v1_c * v2_c)
    except ZeroDivisionError:
        return "Error: Division by zero in relative speed calculation."

    # --- Step 2: Calculate the theoretical total energy ---
    # The total energy of the system is the sum of the individual relativistic energies.
    # E_total = E1 + E2 = gamma1*m1*c^2 + gamma2*m2*c^2
    # We calculate the coefficient of mc^2.

    # Calculate Lorentz factor and energy for astronaut 1
    try:
        gamma1 = 1 / math.sqrt(1 - v1_c**2)
        E1_calc_coeff = gamma1 * m1_coeff
    except ValueError:
        return "Error: Math domain error in Lorentz factor calculation for astronaut 1 (speed >= c)."

    # Calculate Lorentz factor and energy for astronaut 2
    try:
        gamma2 = 1 / math.sqrt(1 - v2_c**2)
        E2_calc_coeff = gamma2 * m2_coeff
    except ValueError:
        return "Error: Math domain error in Lorentz factor calculation for astronaut 2 (speed >= c)."

    # Sum the energies
    E_total_calc_coeff = E1_calc_coeff + E2_calc_coeff

    # --- Step 3: Extract values from the provided answer ---
    # The provided answer is <<<A>>>.
    # Option A is: v_rel = 0.14c , E= 5.96 mc^2
    expected_v_rel_coeff = 0.14
    expected_E_total_coeff = 5.96

    # --- Step 4: Compare calculated values with the answer's values ---
    # The options are given to two decimal places, so we round our calculated values
    # to two decimal places for comparison.
    v_rel_matches = round(v_rel_calc_coeff, 2) == expected_v_rel_coeff
    E_total_matches = round(E_total_calc_coeff, 2) == expected_E_total_coeff

    # --- Step 5: Formulate the final verdict ---
    if v_rel_matches and E_total_matches:
        return "Correct"
    else:
        error_messages = []
        if not v_rel_matches:
            error_messages.append(
                f"The relative speed in the answer is incorrect. "
                f"Expected: {expected_v_rel_coeff}c. "
                f"Calculated: {v_rel_calc_coeff:.4f}c (which rounds to {round(v_rel_calc_coeff, 2)}c)."
            )
        if not E_total_matches:
            error_messages.append(
                f"The total energy in the answer is incorrect. "
                f"Expected: {expected_E_total_coeff}mc^2. "
                f"Calculated: {E_total_calc_coeff:.4f}mc^2 (which rounds to {round(E_total_calc_coeff, 2)}mc^2)."
            )
        return "\n".join(error_messages)

# Execute the check and print the result
result = check_correctness_of_answer()
print(result)