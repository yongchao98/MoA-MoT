import math

def check_answer():
    """
    This function checks the correctness of the given answer for the relativistic physics problem.
    """
    # Given values from the question
    # Velocities are given as fractions of c
    v1_c = 0.6
    v2_c = 0.5
    # Masses are given as multiples of m
    m1_m = 2.0
    m2_m = 3.0

    # --- Step 1: Calculate the theoretical relative speed ---
    # The relativistic velocity subtraction formula for objects moving in the same direction is:
    # v_rel = (v1 - v2) / (1 - (v1 * v2) / c^2)
    # Since velocities are in units of c, the formula simplifies to:
    # v_rel_c = (v1_c - v2_c) / (1 - v1_c * v2_c)
    try:
        calculated_v_rel_c = (v1_c - v2_c) / (1 - v1_c * v2_c)
    except ZeroDivisionError:
        return "Error: Division by zero in relative speed calculation."

    # --- Step 2: Calculate the theoretical total energy ---
    # The total relativistic energy of a particle is E = γmc², where γ is the Lorentz factor.
    # The Lorentz factor γ = 1 / sqrt(1 - (v/c)²).
    # The total energy of the system is E_total = E1 + E2.
    # We will calculate the energy in units of mc².
    
    def calculate_gamma(v_c):
        if abs(v_c) >= 1:
            raise ValueError("Velocity cannot be equal to or greater than the speed of light.")
        return 1 / math.sqrt(1 - v_c**2)

    try:
        # Calculate Lorentz factors for each astronaut
        gamma1 = calculate_gamma(v1_c)
        gamma2 = calculate_gamma(v2_c)

        # Calculate energy for each astronaut in units of mc²
        # E1 = γ1 * m1 * c² = γ1 * (2m) * c² = (γ1 * 2) * mc²
        E1_mc2 = gamma1 * m1_m
        # E2 = γ2 * m2 * c² = γ2 * (3m) * c² = (γ2 * 3) * mc²
        E2_mc2 = gamma2 * m2_m

        # Calculate total system energy in units of mc²
        calculated_E_total_mc2 = E1_mc2 + E2_mc2
    except ValueError as e:
        return f"Error in energy calculation: {e}"


    # --- Step 3: Parse the provided answer and compare ---
    # The provided answer is 'C'.
    # The options are:
    # A) v_rel = 1.1c , E= mc^2
    # B) v_rel=0.14c, E=5mc^2
    # C) v_rel = 0.14c , E= 5.96 mc^2
    # D) v_rel =0.1c , E= 4.96 mc^2
    # So, for option C, the expected values are:
    expected_v_rel_c = 0.14
    expected_E_total_mc2 = 5.96

    # Set a tolerance for floating-point comparison, as the options are rounded.
    # A tolerance of 0.005 is appropriate for values rounded to two decimal places.
    tolerance = 0.005

    # Check if the calculated values match the expected values from option C
    v_rel_correct = abs(calculated_v_rel_c - expected_v_rel_c) < tolerance
    E_total_correct = abs(calculated_E_total_mc2 - expected_E_total_mc2) < tolerance

    if v_rel_correct and E_total_correct:
        return "Correct"
    else:
        error_messages = []
        if not v_rel_correct:
            error_messages.append(
                f"Relative speed does not match. Expected (from option C): {expected_v_rel_c:.2f}c, Calculated: {calculated_v_rel_c:.4f}c."
            )
        if not E_total_correct:
            error_messages.append(
                f"Total energy does not match. Expected (from option C): {expected_E_total_mc2:.2f}mc^2, Calculated: {calculated_E_total_mc2:.4f}mc^2."
            )
        return "Incorrect. " + " ".join(error_messages)

# Run the check and print the result
print(check_answer())