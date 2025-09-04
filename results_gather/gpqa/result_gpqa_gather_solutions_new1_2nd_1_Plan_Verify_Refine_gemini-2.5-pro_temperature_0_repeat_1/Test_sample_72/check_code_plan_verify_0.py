import math

def check_correctness():
    """
    This function checks the correctness of the selected answer for the given physics problem.
    It calculates the correct values for relative speed and total energy and compares them
    to the values in the chosen option.
    """

    # --- Problem Parameters ---
    # Velocities are given as fractions of c.
    v1 = 0.6
    v2 = 0.5
    # Masses are given as factors of m.
    m1_factor = 2.0
    m2_factor = 3.0

    # --- Step 1: Calculate the correct relative speed (v_rel) ---
    # At relativistic speeds, we must use the Einstein velocity subtraction formula:
    # v_rel = (v1 - v2) / (1 - (v1 * v2) / c^2)
    # Since v1 and v2 are in units of c, the formula simplifies to:
    correct_v_rel = (v1 - v2) / (1 - v1 * v2)

    # --- Step 2: Calculate the correct total energy (E_total) ---
    # The total energy of a relativistic particle is E = γmc², where γ is the Lorentz factor.
    # The Lorentz factor γ = 1 / sqrt(1 - v²/c²).
    # The total energy of the system is the sum of the individual energies.
    # We will calculate the numerical factor in front of mc².

    # Lorentz factor for astronaut 1
    gamma1 = 1 / math.sqrt(1 - v1**2)
    # Energy for astronaut 1
    energy1_factor = gamma1 * m1_factor

    # Lorentz factor for astronaut 2
    gamma2 = 1 / math.sqrt(1 - v2**2)
    # Energy for astronaut 2
    energy2_factor = gamma2 * m2_factor

    # Total energy of the system
    correct_energy_factor = energy1_factor + energy2_factor

    # --- Step 3: Define the answer to be checked ---
    # The provided answer is <<<B>>>, which corresponds to:
    # B) v_rel = 0.14c , E= 5.96 mc^2
    answer_v_rel = 0.14
    answer_energy_factor = 5.96

    # --- Step 4: Compare the calculated values with the answer's values ---
    # We use a tolerance because the options are rounded to two decimal places.
    # A tolerance of 0.005 is appropriate for checking values rounded to 2 decimal places.
    tolerance = 0.005

    v_rel_is_correct = math.isclose(correct_v_rel, answer_v_rel, abs_tol=tolerance)
    energy_is_correct = math.isclose(correct_energy_factor, answer_energy_factor, abs_tol=tolerance)

    if v_rel_is_correct and energy_is_correct:
        return "Correct"
    else:
        error_messages = []
        if not v_rel_is_correct:
            error_messages.append(
                f"Constraint violated: Relative Speed. "
                f"The answer's relative speed is {answer_v_rel:.2f}c, but the correct calculation yields {correct_v_rel:.4f}c."
            )
        if not energy_is_correct:
            error_messages.append(
                f"Constraint violated: Total Energy. "
                f"The answer's total energy factor is {answer_energy_factor:.2f}, but the correct calculation yields {correct_energy_factor:.4f}."
            )
        return "\n".join(error_messages)

# Run the check and print the result
result = check_correctness()
print(result)