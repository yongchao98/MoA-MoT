import math

def check_correctness_of_answer():
    """
    This function checks the correctness of the given answer to the relativistic physics problem.
    It calculates the relative speed and total energy from scratch and compares them
    to the values in the selected option.
    """

    # --- Problem Parameters ---
    # Masses are given as multiples of a base mass 'm'
    m1_factor = 2.0
    m2_factor = 3.0
    # Velocities are given as fractions of the speed of light 'c'
    v1_factor = 0.6
    v2_factor = 0.5

    # --- Step 1: Calculate the Relative Speed (v_rel) ---
    # The relativistic velocity subtraction formula for co-linear motion is:
    # v_rel = (v1 - v2) / (1 - (v1 * v2) / c^2)
    # Since we use factors of c, the formula simplifies to:
    v_rel_calculated = (v1_factor - v2_factor) / (1 - v1_factor * v2_factor)

    # --- Step 2: Calculate the Total Energy of the System (E_total) ---
    # The total energy of a particle is E = gamma * m_rest * c^2
    # The Lorentz factor (gamma) is gamma = 1 / sqrt(1 - (v/c)^2)
    # The total energy of the system is E_total = E1 + E2.
    # We will calculate the total energy in units of mc^2.

    # Energy of Astronaut 1
    gamma1 = 1 / math.sqrt(1 - v1_factor**2)
    E1_factor = gamma1 * m1_factor  # This is in units of mc^2

    # Energy of Astronaut 2
    gamma2 = 1 / math.sqrt(1 - v2_factor**2)
    E2_factor = gamma2 * m2_factor  # This is in units of mc^2

    # Total system energy
    E_total_calculated = E1_factor + E2_factor

    # --- Step 3: Verify the selected answer ---
    # The final answer provided is <<<A>>>.
    # The options from the question are:
    # A) v_rel = 0.14c , E= 5.96 mc^2
    # B) v_rel=0.14c, E=5mc^2
    # C) v_rel = 1.1c , E= mc^2
    # D) v_rel =0.1c , E= 4.96 mc^2
    # We need to check if our calculated values match option A.

    expected_v_rel = 0.14
    expected_E = 5.96

    # We use a tolerance for comparing floating-point numbers, as the options are rounded.
    # A tolerance of 0.01 is appropriate since the answers are given to two decimal places.
    v_rel_is_correct = math.isclose(v_rel_calculated, expected_v_rel, abs_tol=0.005)
    E_total_is_correct = math.isclose(E_total_calculated, expected_E, abs_tol=0.005)

    if v_rel_is_correct and E_total_is_correct:
        return "Correct"
    else:
        error_messages = []
        if not v_rel_is_correct:
            error_messages.append(f"The relative speed is incorrect. The calculated value is v_rel ≈ {v_rel_calculated:.4f}c, which rounds to {v_rel_calculated:.2f}c. The answer states it is {expected_v_rel}c.")
        if not E_total_is_correct:
            error_messages.append(f"The total energy is incorrect. The calculated value is E ≈ {E_total_calculated:.4f}mc^2, which rounds to {E_total_calculated:.2f}mc^2. The answer states it is {expected_E}mc^2.")
        return "Incorrect. " + " ".join(error_messages)

# Execute the check and print the result
result = check_correctness_of_answer()
print(result)