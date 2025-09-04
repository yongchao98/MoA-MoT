import math

def check_relativistic_system():
    """
    This function calculates the relative speed and total energy for the given physics problem
    and compares them to the values in the proposed answer (Option D).
    """
    # --- Given parameters from the question ---
    # Velocities are given as fractions of the speed of light, c.
    v1_c = 0.6
    v2_c = 0.5
    # Masses are given as multiples of a base mass, m.
    m1_m = 2.0
    m2_m = 3.0

    # --- Step 1: Calculate the theoretical relative speed ---
    # The relativistic velocity subtraction formula is: v_rel = (v1 - v2) / (1 - (v1*v2)/c^2)
    # We calculate the coefficient that multiplies c.
    try:
        calculated_v_rel_coeff = (v1_c - v2_c) / (1 - v1_c * v2_c)
    except ZeroDivisionError:
        return "Error: Division by zero in relative speed calculation."

    # --- Step 2: Calculate the theoretical total energy ---
    # The total energy is the sum of individual relativistic energies: E_total = E1 + E2
    # The formula for relativistic energy is E = gamma * m_rest * c^2
    # The Lorentz factor is gamma = 1 / sqrt(1 - (v/c)^2)
    # We calculate the coefficient that multiplies mc^2.

    # Energy of the first astronaut
    try:
        gamma1 = 1 / math.sqrt(1 - v1_c**2)
        E1_coeff = gamma1 * m1_m
    except (ValueError, ZeroDivisionError) as e:
        return f"Error calculating energy for astronaut 1: {e}"

    # Energy of the second astronaut
    try:
        gamma2 = 1 / math.sqrt(1 - v2_c**2)
        E2_coeff = gamma2 * m2_m
    except (ValueError, ZeroDivisionError) as e:
        return f"Error calculating energy for astronaut 2: {e}"

    calculated_E_total_coeff = E1_coeff + E2_coeff

    # --- Step 3: Define the proposed answer's values ---
    # The final answer provided is <<<D>>>, which corresponds to:
    # D) v_rel = 0.14c , E= 5.96 mc^2
    proposed_v_rel_coeff = 0.14
    proposed_E_total_coeff = 5.96

    # --- Step 4: Compare the calculated values with the proposed answer ---
    # The options are given to two decimal places, so we compare our more precise
    # calculations by rounding them to the same precision.
    v_rel_matches = round(calculated_v_rel_coeff, 2) == proposed_v_rel_coeff
    E_total_matches = round(calculated_E_total_coeff, 2) == proposed_E_total_coeff

    if v_rel_matches and E_total_matches:
        return "Correct"
    else:
        error_message = "Incorrect. The values in the proposed answer do not match the calculated results.\n"
        if not v_rel_matches:
            error_message += (f"Constraint failure (Relative Speed): "
                              f"Calculated v_rel is approximately {calculated_v_rel_coeff:.4f}c, which rounds to {round(calculated_v_rel_coeff, 2)}c. "
                              f"The proposed answer has v_rel = {proposed_v_rel_coeff}c.\n")
        if not E_total_matches:
            error_message += (f"Constraint failure (Total Energy): "
                              f"Calculated E_total is approximately {calculated_E_total_coeff:.4f}mc^2, which rounds to {round(calculated_E_total_coeff, 2)}mc^2. "
                              f"The proposed answer has E_total = {proposed_E_total_coeff}mc^2.\n")
        return error_message.strip()

# Execute the check and print the result
result = check_relativistic_system()
print(result)