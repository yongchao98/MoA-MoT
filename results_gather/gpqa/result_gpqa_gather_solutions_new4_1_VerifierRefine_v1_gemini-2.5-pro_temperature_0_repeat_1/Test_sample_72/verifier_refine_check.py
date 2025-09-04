import math

def check_correctness_of_answer():
    """
    This function checks the correctness of the provided answer to the relativistic physics problem.
    It recalculates the relative speed and total energy based on the problem statement and compares
    them to the values given in the selected answer (Option A).
    """

    # --- Problem Parameters from the Question ---
    # Astronaut 1
    m1_factor = 2.0  # mass in units of 'm'
    v1_factor = 0.6  # velocity as a fraction of 'c'

    # Astronaut 2
    m2_factor = 3.0  # mass in units of 'm'
    v2_factor = 0.5  # velocity as a fraction of 'c'

    # --- Values from the Proposed Answer (Option A) ---
    # The final answer provided is <<<A>>>, which corresponds to:
    # v_rel = 0.14c , E = 5.96 mc^2
    v_rel_answer = 0.14
    E_total_answer_factor = 5.96

    # --- Step 1: Recalculate the Relative Speed (v_rel) ---
    # The relativistic velocity subtraction formula for co-linear motion is:
    # v_rel = (v1 - v2) / (1 - (v1 * v2) / c^2)
    # Since velocities are given as factors of c, the formula simplifies.
    try:
        v_rel_calculated = (v1_factor - v2_factor) / (1 - v1_factor * v2_factor)
    except ZeroDivisionError:
        return "Calculation Error: Division by zero occurred while calculating relative speed."

    # --- Step 2: Recalculate the Total Energy (E_total) ---
    # The total energy of the system is the sum of individual relativistic energies:
    # E_total = E1 + E2 = (gamma1*m1 + gamma2*m2) * c^2
    # where the Lorentz factor gamma = 1 / sqrt(1 - (v/c)^2)
    try:
        # Calculate energy for astronaut 1
        gamma1 = 1 / math.sqrt(1 - v1_factor**2)
        E1_factor = gamma1 * m1_factor

        # Calculate energy for astronaut 2
        gamma2 = 1 / math.sqrt(1 - v2_factor**2)
        E2_factor = gamma2 * m2_factor

        # Calculate the total energy factor (the coefficient of mc^2)
        E_total_calculated_factor = E1_factor + E2_factor
    except ValueError:
        return "Calculation Error: Math domain error, likely due to a speed >= c."

    # --- Step 3: Compare Calculated Values with the Answer's Values ---
    # The options are given to two decimal places, so we round our calculated
    # values to the same precision for a fair comparison.
    is_v_rel_correct = (round(v_rel_calculated, 2) == v_rel_answer)
    is_E_total_correct = (round(E_total_calculated_factor, 2) == E_total_answer_factor)

    if is_v_rel_correct and is_E_total_correct:
        return "Correct"
    else:
        error_messages = []
        if not is_v_rel_correct:
            error_messages.append(
                f"The relative speed is incorrect. The calculated value is {v_rel_calculated:.4f}c, which rounds to {round(v_rel_calculated, 2)}c, but the answer states {v_rel_answer}c."
            )
        if not is_E_total_correct:
            error_messages.append(
                f"The total energy is incorrect. The calculated value is {E_total_calculated_factor:.4f}mc^2, which rounds to {round(E_total_calculated_factor, 2)}mc^2, but the answer states {E_total_answer_factor}mc^2."
            )
        return "Incorrect. " + " ".join(error_messages)

# Execute the check and print the result
result = check_correctness_of_answer()
print(result)