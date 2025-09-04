import math

def check_relativity_problem():
    """
    This function checks the correctness of the answer to the relativistic physics problem.
    It calculates the relative speed and total energy from scratch and compares them
    to the values given in the selected answer option.
    """
    # --- Problem Parameters ---
    # Velocities are given as fractions of the speed of light, c.
    v1 = 0.6  # Velocity of astronaut 1
    v2 = 0.5  # Velocity of astronaut 2
    # Masses are given as multiples of a base mass, m.
    m1_factor = 2  # Mass of astronaut 1 is 2m
    m2_factor = 3  # Mass of astronaut 2 is 3m

    # --- Candidate Answer to Check ---
    # The provided answer is <<<D>>>, which corresponds to:
    # D) v_rel = 0.14c , E= 5.96 mc^2
    # We extract these numerical values for comparison.
    expected_v_rel = 0.14  # in units of c
    expected_E_total = 5.96  # in units of mc^2

    # --- Calculation Step 1: Relative Speed (v_rel) ---
    # The relativistic velocity subtraction formula for objects moving in the same direction is:
    # v_rel = (v1 - v2) / (1 - (v1 * v2) / c^2)
    # Since v1 and v2 are already in units of c, the formula simplifies to:
    try:
        calculated_v_rel = (v1 - v2) / (1 - v1 * v2)
    except ZeroDivisionError:
        return "Error: Division by zero during relative speed calculation."

    # --- Calculation Step 2: Total Energy (E_total) ---
    # The total energy of a relativistic particle is E = γmc², where γ is the Lorentz factor.
    # The Lorentz factor is γ = 1 / sqrt(1 - v²/c²).
    # The total energy of the system is the sum of the individual energies: E_total = E1 + E2.
    # We calculate the total energy in units of mc².
    try:
        # Energy of astronaut 1
        gamma1 = 1 / math.sqrt(1 - v1**2)
        E1 = gamma1 * m1_factor  # Energy in units of mc^2

        # Energy of astronaut 2
        gamma2 = 1 / math.sqrt(1 - v2**2)
        E2 = gamma2 * m2_factor  # Energy in units of mc^2

        # Total energy of the system
        calculated_E_total = E1 + E2
    except ValueError:
        return "Error: Math domain error during energy calculation (e.g., speed > c)."

    # --- Verification Step ---
    # The options in the question are rounded to two decimal places.
    # We will compare our more precise calculated values by rounding them to two decimal places
    # and checking if they match the expected answer.
    
    # Check if the relative speed matches
    v_rel_matches = math.isclose(round(calculated_v_rel, 2), expected_v_rel, rel_tol=1e-9)
    if not v_rel_matches:
        return (f"Incorrect. The relative speed constraint is not satisfied. "
                f"The calculated relative speed is {calculated_v_rel:.4f}c, which rounds to {round(calculated_v_rel, 2)}c. "
                f"The answer D states v_rel = {expected_v_rel}c.")

    # Check if the total energy matches
    E_total_matches = math.isclose(round(calculated_E_total, 2), expected_E_total, rel_tol=1e-9)
    if not E_total_matches:
        return (f"Incorrect. The total energy constraint is not satisfied. "
                f"The calculated total energy is {calculated_E_total:.4f}mc^2, which rounds to {round(calculated_E_total, 2)}mc^2. "
                f"The answer D states E = {expected_E_total}mc^2.")

    # If both values match the expected answer
    return "Correct"

# Run the check and print the result
result = check_relativity_problem()
print(result)