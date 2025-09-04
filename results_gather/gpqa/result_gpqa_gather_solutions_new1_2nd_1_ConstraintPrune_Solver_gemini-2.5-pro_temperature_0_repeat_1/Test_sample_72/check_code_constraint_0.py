import math

def check_correctness_of_answer():
    """
    This function checks the correctness of the final answer for the given physics problem.
    It calculates the relative speed and total energy based on the principles of special relativity
    and compares them to the values in the selected option.
    """

    # --- Problem Parameters ---
    # Astronaut 1: mass = 2m, velocity = 0.6c
    # Astronaut 2: mass = 3m, velocity = 0.5c
    m1_factor = 2.0
    v1_factor = 0.6
    m2_factor = 3.0
    v2_factor = 0.5

    # --- Proposed Answer ---
    # The final answer given is 'D', which corresponds to:
    # v_rel = 0.14c , E = 5.96 mc^2
    # We extract the numerical factors to check against.
    expected_v_rel_factor = 0.14
    expected_E_factor = 5.96

    # --- Calculation from First Principles ---

    # 1. Calculate the relative speed (v_rel)
    # The relativistic velocity subtraction formula for objects moving in the same direction is:
    # v_rel = (v1 - v2) / (1 - (v1 * v2) / c^2)
    # When working with factors of c, the formula simplifies to:
    calculated_v_rel_factor = (v1_factor - v2_factor) / (1 - v1_factor * v2_factor)

    # 2. Calculate the total energy (E_total)
    # The total energy of the system is the sum of the individual total relativistic energies:
    # E_total = E1 + E2 = γ1*m1*c^2 + γ2*m2*c^2
    # We calculate the numerical factor in front of mc^2.

    def lorentz_factor(v_factor):
        """Calculates the Lorentz factor γ for a given velocity factor v/c."""
        if abs(v_factor) >= 1:
            return float('inf')  # Unphysical speed
        return 1 / math.sqrt(1 - v_factor**2)

    # Calculate Lorentz factors for each astronaut
    gamma1 = lorentz_factor(v1_factor)
    gamma2 = lorentz_factor(v2_factor)

    # Calculate the energy factor (in units of mc^2) for each astronaut
    # E1_factor = γ1 * (mass_factor_1)
    # E2_factor = γ2 * (mass_factor_2)
    E1_factor = gamma1 * m1_factor
    E2_factor = gamma2 * m2_factor

    # The total energy factor is the sum of the individual energy factors
    calculated_E_factor = E1_factor + E2_factor

    # --- Verification ---
    # The options in the question are rounded to two decimal places.
    # We must compare our calculated values, rounded to the same precision, with the expected values.

    # Check if the rounded calculated relative speed matches the expected value.
    v_rel_matches = math.isclose(round(calculated_v_rel_factor, 2), expected_v_rel_factor, rel_tol=1e-9)

    # Check if the rounded calculated total energy matches the expected value.
    E_matches = math.isclose(round(calculated_E_factor, 2), expected_E_factor, rel_tol=1e-9)

    # --- Return Result ---
    if v_rel_matches and E_matches:
        return "Correct"
    else:
        error_messages = []
        if not v_rel_matches:
            error_messages.append(
                f"The relative speed is incorrect. The answer states v_rel = {expected_v_rel_factor}c, "
                f"but the calculated value is approximately {calculated_v_rel_factor:.4f}c, which rounds to {round(calculated_v_rel_factor, 2)}c."
            )
        if not E_matches:
            error_messages.append(
                f"The total energy is incorrect. The answer states E = {expected_E_factor}mc^2, "
                f"but the calculated value is approximately {calculated_E_factor:.4f}mc^2, which rounds to {round(calculated_E_factor, 2)}mc^2."
            )
        return "Incorrect. " + " ".join(error_messages)

# Execute the check and print the result
result = check_correctness_of_answer()
print(result)