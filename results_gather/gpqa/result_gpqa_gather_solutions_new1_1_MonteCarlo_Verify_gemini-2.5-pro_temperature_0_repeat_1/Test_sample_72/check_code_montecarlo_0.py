import math

def check_correctness_of_answer():
    """
    This function checks the correctness of the provided answer to the physics problem.
    It calculates the theoretical values for relative speed and total energy and compares them
    to the values given in the selected answer option 'A'.
    """

    # --- Problem Parameters ---
    # Astronaut 1
    m1_factor = 2.0  # Mass in units of 'm'
    v1_factor = 0.6  # Velocity in units of 'c'

    # Astronaut 2
    m2_factor = 3.0  # Mass in units of 'm'
    v2_factor = 0.5  # Velocity in units of 'c'

    # --- Answer to be Checked (Option A) ---
    # The final answer provided is 'A', which corresponds to:
    # v_rel = 0.14c, E = 5.96 mc^2
    expected_v_rel_factor = 0.14
    expected_E_factor = 5.96

    # --- Calculation ---

    # 1. Calculate the theoretical relative speed (v_rel)
    # The relativistic velocity subtraction formula for objects moving in the same direction is:
    # v_rel = (v1 - v2) / (1 - (v1 * v2) / c^2)
    # We can calculate the factor by which c is multiplied:
    try:
        calculated_v_rel_factor = (v1_factor - v2_factor) / (1 - v1_factor * v2_factor)
    except ZeroDivisionError:
        return "Calculation Error: Division by zero when calculating relative speed."

    # 2. Calculate the theoretical total energy (E_total)
    # The total energy of the system is the sum of the individual relativistic energies.
    # E_total = E1 + E2
    # The formula for a single particle's energy is E = gamma * m * c^2,
    # where gamma = 1 / sqrt(1 - (v/c)^2).
    try:
        # Calculate gamma and energy for astronaut 1
        gamma1 = 1 / math.sqrt(1 - v1_factor**2)
        E1_factor = gamma1 * m1_factor  # Energy in units of mc^2

        # Calculate gamma and energy for astronaut 2
        gamma2 = 1 / math.sqrt(1 - v2_factor**2)
        E2_factor = gamma2 * m2_factor  # Energy in units of mc^2

        # Calculate total energy factor
        calculated_E_factor = E1_factor + E2_factor
    except ValueError:
        return "Calculation Error: Math domain error (e.g., sqrt of a negative number) during energy calculation. A velocity might be >= c."
    except ZeroDivisionError:
        return "Calculation Error: Division by zero during energy calculation. A velocity might be equal to c."


    # --- Verification ---
    # Compare the calculated values with the expected values from the answer,
    # using a tolerance to account for rounding in the options.
    tolerance = 0.01

    v_rel_is_correct = abs(calculated_v_rel_factor - expected_v_rel_factor) < tolerance
    E_is_correct = abs(calculated_E_factor - expected_E_factor) < tolerance

    if v_rel_is_correct and E_is_correct:
        return "Correct"
    else:
        error_messages = []
        if not v_rel_is_correct:
            error_messages.append(
                f"Relative speed constraint is not satisfied. "
                f"Expected approximately {expected_v_rel_factor}c, but the calculated value is {calculated_v_rel_factor:.4f}c."
            )
        if not E_is_correct:
            error_messages.append(
                f"Total energy constraint is not satisfied. "
                f"Expected approximately {expected_E_factor}mc^2, but the calculated value is {calculated_E_factor:.4f}mc^2."
            )
        return "\n".join(error_messages)

# Run the check and print the result
print(check_correctness_of_answer())