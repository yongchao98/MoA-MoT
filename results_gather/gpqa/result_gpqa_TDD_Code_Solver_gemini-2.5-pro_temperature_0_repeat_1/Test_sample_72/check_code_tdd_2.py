import math

def check_relativistic_kinematics():
    """
    This function checks the correctness of the provided answer for a relativistic physics problem.
    It calculates the relative speed and total energy of a system of two astronauts and compares
    the results with the given answer option.
    """
    # --- Problem Parameters ---
    # Masses are in units of 'm', speeds are in units of 'c'
    m1_factor = 2.0
    v1 = 0.6
    
    m2_factor = 3.0
    v2 = 0.5

    # --- Answer to Check (Option D) ---
    # Relative speed is in units of 'c', Energy is in units of 'mc^2'
    answer_v_rel = 0.14
    answer_E_total = 5.96

    # --- Calculation ---

    # 1. Calculate Relative Speed
    # The relativistic velocity subtraction formula is v_rel = (v_b - v_a) / (1 - (v_a * v_b) / c^2).
    # Since speeds are given as fractions of c, the formula simplifies.
    # The question asks for speed, which is the magnitude of the velocity.
    try:
        v_rel_velocity = (v2 - v1) / (1 - v1 * v2)
        calculated_v_rel_speed = abs(v_rel_velocity)
    except Exception as e:
        return f"An error occurred during relative speed calculation: {e}"

    # 2. Calculate Total Energy
    # The total energy of a particle is E = gamma * m_rest * c^2.
    # The Lorentz factor (gamma) is 1 / sqrt(1 - v^2/c^2).
    # The total energy of the system is the sum of the individual energies.
    try:
        # Astronaut 1
        gamma1 = 1 / math.sqrt(1 - v1**2)
        E1_factor = gamma1 * m1_factor  # This is the coefficient of mc^2

        # Astronaut 2
        gamma2 = 1 / math.sqrt(1 - v2**2)
        E2_factor = gamma2 * m2_factor  # This is the coefficient of mc^2

        # Total Energy
        calculated_E_total_factor = E1_factor + E2_factor
    except Exception as e:
        return f"An error occurred during energy calculation: {e}"

    # --- Verification ---
    # We compare the calculated values with the answer's values.
    # We use a tolerance because the answer is rounded to two decimal places.
    # A relative tolerance of 1% is appropriate for checking values rounded to 2 significant figures.
    
    # Check relative speed
    if not math.isclose(calculated_v_rel_speed, answer_v_rel, rel_tol=0.01):
        return (f"Incorrect: The relative speed is calculated to be {calculated_v_rel_speed:.4f}c, "
                f"which rounds to {calculated_v_rel_speed:.2f}c. The answer provided is {answer_v_rel}c, "
                f"which does not match the calculated value.")

    # Check total energy
    if not math.isclose(calculated_E_total_factor, answer_E_total, rel_tol=0.001):
        return (f"Incorrect: The total energy is calculated to be {calculated_E_total_factor:.4f}mc^2, "
                f"which rounds to {calculated_E_total_factor:.2f}mc^2. The answer provided is {answer_E_total}mc^2, "
                f"which does not match the calculated value.")

    # If both checks pass, the answer is correct.
    return "Correct"

# Execute the check and print the result
result = check_relativistic_kinematics()
print(result)