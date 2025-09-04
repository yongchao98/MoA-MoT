import math

def check_correctness_of_answer():
    """
    This function checks the correctness of the final answer for the given physics problem.
    It calculates the relative speed and total energy from scratch and compares them
    to the values in the selected option 'D'.
    """

    # --- Problem Parameters ---
    # Astronaut 1
    m1_m = 2.0  # Mass in units of 'm'
    v1_c = 0.6  # Velocity in units of 'c'

    # Astronaut 2
    m2_m = 3.0  # Mass in units of 'm'
    v2_c = 0.5  # Velocity in units of 'c'

    # --- Step 1: Calculate the Relative Speed (v_rel) ---
    # The relativistic velocity subtraction formula for objects moving in the same direction is:
    # v_rel = (v1 - v2) / (1 - (v1 * v2) / c^2)
    # In terms of c, this is: v_rel_c = (v1_c - v2_c) / (1 - v1_c * v2_c)
    try:
        calculated_v_rel_c = (v1_c - v2_c) / (1 - v1_c * v2_c)
    except ZeroDivisionError:
        return "Calculation Error: Division by zero in relative speed calculation."

    # --- Step 2: Calculate the Total Energy (E_total) ---
    # The total energy of the system is the sum of the individual relativistic energies.
    # E_total = E1 + E2 = (gamma1 * m1 * c^2) + (gamma2 * m2 * c^2)
    # The Lorentz factor (gamma) is: 1 / sqrt(1 - (v/c)^2)
    # We will calculate the total energy in units of mc^2.

    try:
        # Energy of Astronaut 1
        gamma1 = 1 / math.sqrt(1 - v1_c**2)
        E1_mc2 = gamma1 * m1_m

        # Energy of Astronaut 2
        gamma2 = 1 / math.sqrt(1 - v2_c**2)
        E2_mc2 = gamma2 * m2_m

        # Total Energy of the system
        calculated_E_total_mc2 = E1_mc2 + E2_mc2
    except (ValueError, ZeroDivisionError) as e:
        return f"Calculation Error: {e} in energy calculation."

    # --- Step 3: Verify the Final Answer (Option D) ---
    # The provided final answer is D: v_rel = 0.14c, E = 5.96 mc^2
    target_v_rel_c = 0.14
    target_E_total_mc2 = 5.96

    # Round the calculated values to two decimal places to match the format of the options.
    v_rel_rounded = round(calculated_v_rel_c, 2)
    E_total_rounded = round(calculated_E_total_mc2, 2)

    # Check if the rounded calculated values match the values in option D.
    v_rel_match = (v_rel_rounded == target_v_rel_c)
    E_total_match = (E_total_rounded == target_E_total_mc2)

    if v_rel_match and E_total_match:
        return "Correct"
    else:
        error_messages = []
        if not v_rel_match:
            error_messages.append(
                f"The relative speed in option D is incorrect. "
                f"Expected v_rel ≈ {v_rel_rounded}c (from calculation {calculated_v_rel_c:.4f}c), "
                f"but option D has v_rel = {target_v_rel_c}c."
            )
        if not E_total_match:
            error_messages.append(
                f"The total energy in option D is incorrect. "
                f"Expected E ≈ {E_total_rounded}mc^2 (from calculation {calculated_E_total_mc2:.4f}mc^2), "
                f"but option D has E = {target_E_total_mc2}mc^2."
            )
        return " ".join(error_messages)

# Execute the check and print the result
result = check_correctness_of_answer()
print(result)