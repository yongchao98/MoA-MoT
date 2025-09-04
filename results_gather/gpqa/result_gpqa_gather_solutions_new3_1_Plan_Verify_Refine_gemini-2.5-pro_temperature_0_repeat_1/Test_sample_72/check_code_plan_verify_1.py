import math

def check_correctness():
    """
    This function checks the correctness of the provided answer to the physics problem.
    It calculates the relative speed and total energy based on the principles of special relativity
    and compares them to the values in the selected option 'D'.
    """

    # --- Problem Parameters ---
    # Velocities are given as coefficients of c (speed of light)
    v1_coeff = 0.6
    v2_coeff = 0.5
    # Masses are given as coefficients of m
    m1_coeff = 2.0
    m2_coeff = 3.0

    # --- 1. Calculate the theoretical relative speed ---
    # The relativistic velocity subtraction formula for co-linear motion is:
    # v_rel = (v1 - v2) / (1 - (v1 * v2) / c^2)
    # We calculate the coefficient of c.
    try:
        calc_v_rel_coeff = (v1_coeff - v2_coeff) / (1 - v1_coeff * v2_coeff)
    except Exception as e:
        return f"An error occurred during relative speed calculation: {e}"

    # --- 2. Calculate the theoretical total energy ---
    # The total energy of the system is the sum of the individual relativistic energies.
    # E_total = E1 + E2, where E = gamma * m * c^2
    # gamma = 1 / sqrt(1 - (v/c)^2)
    # We calculate the coefficient of mc^2.
    try:
        # Energy of astronaut 1
        gamma1 = 1 / math.sqrt(1 - v1_coeff**2)
        E1_coeff = gamma1 * m1_coeff

        # Energy of astronaut 2
        gamma2 = 1 / math.sqrt(1 - v2_coeff**2)
        E2_coeff = gamma2 * m2_coeff

        # Total energy coefficient
        calc_E_total_coeff = E1_coeff + E2_coeff
    except Exception as e:
        return f"An error occurred during total energy calculation: {e}"

    # --- 3. Extract values from the selected answer ---
    # The final answer given is 'D', which corresponds to:
    # v_rel = 0.14c , E = 5.96 mc^2
    expected_v_rel_coeff = 0.14
    expected_E_total_coeff = 5.96

    # --- 4. Compare calculated values with the answer's values ---
    # A tolerance is needed because the options are rounded to two decimal places.
    # A tolerance of 0.005 is appropriate for values rounded to 2 decimal places.
    tolerance = 0.005

    # Check relative speed
    if not math.isclose(calc_v_rel_coeff, expected_v_rel_coeff, rel_tol=0, abs_tol=tolerance):
        return (f"The relative speed is incorrect. The answer states v_rel = {expected_v_rel_coeff}c, "
                f"but the calculated value is approximately {calc_v_rel_coeff:.4f}c.")

    # Check total energy
    if not math.isclose(calc_E_total_coeff, expected_E_total_coeff, rel_tol=0, abs_tol=tolerance):
        return (f"The total energy is incorrect. The answer states E = {expected_E_total_coeff}mc^2, "
                f"but the calculated value is approximately {calc_E_total_coeff:.4f}mc^2.")

    # If both checks pass, the answer is correct.
    return "Correct"

# Execute the check and print the result
result = check_correctness()
print(result)