import math

def check_correctness_of_answer():
    """
    This function calculates the relative speed and total energy for the given physics problem
    and compares the results with the values in the selected answer (Option D).
    """
    
    # --- Problem Parameters from the question ---
    # Astronaut 1
    m1_factor = 2.0  # Mass in units of 'm'
    v1_factor = 0.6  # Velocity as a fraction of 'c'
    
    # Astronaut 2
    m2_factor = 3.0  # Mass in units of 'm'
    v2_factor = 0.5  # Velocity as a fraction of 'c'

    # --- Values from the selected answer to be checked (Option D) ---
    # D) v_rel = 0.14c , E= 5.96 mc^2
    expected_v_rel_factor = 0.14
    expected_E_total_factor = 5.96

    # --- Step 1: Calculate the theoretical relative speed ---
    # The relativistic velocity subtraction formula for objects moving in the same direction is:
    # v_rel = (v1 - v2) / (1 - (v1 * v2) / c^2)
    # Since velocities are given as fractions of c, the formula simplifies to:
    v_rel_numerator = v1_factor - v2_factor
    v_rel_denominator = 1 - (v1_factor * v2_factor)
    calculated_v_rel_factor = v_rel_numerator / v_rel_denominator

    # --- Step 2: Calculate the theoretical total energy ---
    # The total energy of the system is the sum of the individual relativistic energies.
    # E_total = E1 + E2
    # The formula for relativistic energy is E = gamma * m_rest * c^2,
    # where the Lorentz factor gamma = 1 / sqrt(1 - (v/c)^2).

    # Calculate energy for Astronaut 1
    gamma1 = 1 / math.sqrt(1 - v1_factor**2)
    E1_factor = gamma1 * m1_factor  # Energy in units of mc^2

    # Calculate energy for Astronaut 2
    gamma2 = 1 / math.sqrt(1 - v2_factor**2)
    E2_factor = gamma2 * m2_factor  # Energy in units of mc^2

    # Calculate total energy of the system
    calculated_E_total_factor = E1_factor + E2_factor

    # --- Step 3: Verify the answer ---
    # We check if the calculated values, when rounded to two decimal places (as in the options),
    # match the values from the selected answer.
    
    v_rel_matches = round(calculated_v_rel_factor, 2) == expected_v_rel_factor
    E_total_matches = round(calculated_E_total_factor, 2) == expected_E_total_factor

    if v_rel_matches and E_total_matches:
        return "Correct"
    else:
        error_messages = []
        if not v_rel_matches:
            error_messages.append(
                f"The relative speed constraint is not satisfied. "
                f"The answer states v_rel = {expected_v_rel_factor}c, but the calculated value is approximately {calculated_v_rel_factor:.4f}c, "
                f"which rounds to {round(calculated_v_rel_factor, 2)}c."
            )
        if not E_total_matches:
            # Check for a common mistake (using only rest mass)
            rest_energy_total = m1_factor + m2_factor
            if math.isclose(expected_E_total_factor, rest_energy_total):
                 error_messages.append(
                    f"The total energy constraint is not satisfied. "
                    f"The answer states E = {expected_E_total_factor}mc^2, which incorrectly represents only the total rest mass energy. "
                    f"The correct total energy, including relativistic kinetic energy, is approximately {calculated_E_total_factor:.4f}mc^2."
                )
            else:
                error_messages.append(
                    f"The total energy constraint is not satisfied. "
                    f"The answer states E = {expected_E_total_factor}mc^2, but the calculated value is approximately {calculated_E_total_factor:.4f}mc^2, "
                    f"which rounds to {round(calculated_E_total_factor, 2)}mc^2."
                )
        
        return "Incorrect. " + " ".join(error_messages)

# Run the check and print the result.
result = check_correctness_of_answer()
print(result)