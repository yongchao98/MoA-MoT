import math

def check_correctness():
    """
    This function checks the correctness of the selected answer for the given physics problem.
    It calculates the theoretical values for relative speed and total energy and compares them
    with the values from the chosen option 'D'.
    """
    
    # --- Problem Parameters ---
    # Velocities are given as fractions of the speed of light, c.
    # Masses are given as multiples of a base mass, m.
    v1_c = 0.6  # Velocity of astronaut 1 in units of c
    v2_c = 0.5  # Velocity of astronaut 2 in units of c
    m1_m = 2.0  # Mass of astronaut 1 in units of m
    m2_m = 3.0  # Mass of astronaut 2 in units of m

    # --- 1. Calculate the Theoretical Relative Speed ---
    # The relativistic velocity subtraction formula for co-linear motion is:
    # v_rel = (v1 - v2) / (1 - (v1 * v2) / c^2)
    # Since we are working with velocities as fractions of c, the formula simplifies to:
    # v_rel_c = (v1_c - v2_c) / (1 - v1_c * v2_c)
    
    calculated_v_rel_c = (v1_c - v2_c) / (1 - v1_c * v2_c)

    # --- 2. Calculate the Theoretical Total Energy ---
    # The total energy of the system is the sum of the individual total relativistic energies.
    # E_total = E1 + E2 = (gamma1 * m1 * c^2) + (gamma2 * m2 * c^2)
    # We will calculate the total energy in units of mc^2.
    # The Lorentz factor (gamma) is given by: gamma = 1 / sqrt(1 - (v/c)^2)

    # Energy of astronaut 1
    gamma1 = 1 / math.sqrt(1 - v1_c**2)
    E1_mc2 = gamma1 * m1_m

    # Energy of astronaut 2
    gamma2 = 1 / math.sqrt(1 - v2_c**2)
    E2_mc2 = gamma2 * m2_m

    # Total energy of the system
    calculated_E_total_mc2 = E1_mc2 + E2_mc2

    # --- 3. Extract Values from the Provided Answer ---
    # The final answer provided is 'D'.
    # Option D is: v_rel = 0.14c, E = 5.96 mc^2
    answer_v_rel_c = 0.14
    answer_E_total_mc2 = 5.96

    # --- 4. Compare Calculated Values with the Answer ---
    # The values in the options are rounded to two decimal places.
    # We use math.isclose() to account for potential floating-point inaccuracies and rounding.
    # A relative tolerance of 1% for velocity and 0.1% for energy is appropriate.
    
    v_rel_is_correct = math.isclose(calculated_v_rel_c, answer_v_rel_c, rel_tol=0.01)
    E_total_is_correct = math.isclose(calculated_E_total_mc2, answer_E_total_mc2, rel_tol=0.001)

    if not v_rel_is_correct:
        return (f"Incorrect: The relative speed in the answer is wrong. "
                f"The calculated relative speed is approximately {calculated_v_rel_c:.4f}c, "
                f"which rounds to {round(calculated_v_rel_c, 2)}c. The answer provides {answer_v_rel_c}c.")

    if not E_total_is_correct:
        return (f"Incorrect: The total energy in the answer is wrong. "
                f"The calculated total energy is approximately {calculated_E_total_mc2:.4f}mc^2, "
                f"which rounds to {round(calculated_E_total_mc2, 2)}mc^2. The answer provides {answer_E_total_mc2}mc^2.")

    return "Correct"

# Run the check and print the result.
result = check_correctness()
print(result)