import math

def check_relativistic_calculation():
    """
    This function checks the correctness of the answer to the given physics problem.
    It calculates the relative speed and total energy of the system and compares
    them to the values in the proposed answer (Option C).
    """
    # --- Problem Parameters ---
    # We can work with factors of m and c, effectively setting m=1 and c=1.
    m1_factor = 2.0  # Mass of astronaut 1 is 2m
    v1_factor = 0.6  # Speed of astronaut 1 is 0.6c
    m2_factor = 3.0  # Mass of astronaut 2 is 3m
    v2_factor = 0.5  # Speed of astronaut 2 is 0.5c

    # --- Answer to Check (Option C) ---
    # The values in the options are rounded to two decimal places.
    expected_v_rel_factor = 0.14
    expected_E_total_factor = 5.96

    # --- 1. Calculate Relative Speed ---
    # The relativistic velocity subtraction formula for collinear motion is:
    # v_rel = (v_a - v_b) / (1 - (v_a * v_b) / c^2)
    # The question asks for speed, which is the magnitude.
    try:
        calc_v_rel_factor = abs(v1_factor - v2_factor) / (1 - (v1_factor * v2_factor))
    except ZeroDivisionError:
        return "Error: Division by zero in relative speed calculation. This would imply v1*v2 = c^2."

    # --- 2. Calculate Total Energy ---
    # The total relativistic energy of a particle is E = gamma * m_rest * c^2,
    # where gamma = 1 / sqrt(1 - v^2/c^2).
    # The total energy of the system is the sum of the individual energies.

    # Calculate gamma and energy for astronaut 1
    try:
        gamma1 = 1 / math.sqrt(1 - v1_factor**2)
        E1_factor = gamma1 * m1_factor  # Energy in units of mc^2
    except (ValueError, ZeroDivisionError) as e:
        return f"Error calculating energy for astronaut 1: {e}. Speed might be >= c."

    # Calculate gamma and energy for astronaut 2
    try:
        gamma2 = 1 / math.sqrt(1 - v2_factor**2)
        E2_factor = gamma2 * m2_factor  # Energy in units of mc^2
    except (ValueError, ZeroDivisionError) as e:
        return f"Error calculating energy for astronaut 2: {e}. Speed might be >= c."

    # Calculate total energy of the system
    calc_E_total_factor = E1_factor + E2_factor

    # --- 3. Verify the Answer ---
    # We check if the calculated values, when rounded to two decimal places,
    # match the values given in the answer option.
    v_rel_correct = round(calc_v_rel_factor, 2) == expected_v_rel_factor
    E_total_correct = round(calc_E_total_factor, 2) == expected_E_total_factor

    if v_rel_correct and E_total_correct:
        return "Correct"
    else:
        error_messages = []
        if not v_rel_correct:
            error_messages.append(
                f"Constraint failed: Relative speed is incorrect. "
                f"The calculated relative speed is {calc_v_rel_factor:.4f}c, which rounds to {round(calc_v_rel_factor, 2)}c. "
                f"The answer states it is {expected_v_rel_factor}c."
            )
        if not E_total_correct:
            error_messages.append(
                f"Constraint failed: Total energy is incorrect. "
                f"The calculated total energy is {calc_E_total_factor:.4f}mc^2, which rounds to {round(calc_E_total_factor, 2)}mc^2. "
                f"The answer states it is {expected_E_total_factor}mc^2."
            )
        return "\n".join(error_messages)

# Run the check and print the result
result = check_relativistic_calculation()
print(result)