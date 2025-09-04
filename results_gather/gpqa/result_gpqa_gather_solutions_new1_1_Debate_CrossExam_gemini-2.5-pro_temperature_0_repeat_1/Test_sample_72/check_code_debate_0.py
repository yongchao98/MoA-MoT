import math

def check_correctness_of_answer():
    """
    This function checks the correctness of the provided answer to the physics problem.
    It calculates the relative speed and total energy based on the principles of special relativity
    and compares them to the values given in the selected option.
    """
    
    # --- Problem Parameters ---
    # We can set c=1 and m=1 for simplicity, as the final answer is in terms of c and mc^2.
    m1_factor = 2.0
    v1 = 0.6  # in units of c
    m2_factor = 3.0
    v2 = 0.5  # in units of c

    # --- The Answer to Check ---
    # The final answer provided is <<<C>>>.
    # The options listed in the final answer's analysis are:
    # A) v_rel =0.1c , E= 4.96 mc^2
    # B) v_rel = 1.1c , E= mc^2
    # C) v_rel = 0.14c , E= 5.96 mc^2
    # D) v_rel=0.14c, E=5mc^2
    # We will check if the calculated values match option C.
    
    answer_v_rel = 0.14  # from option C
    answer_E_total = 5.96 # from option C

    # --- Step 1: Calculate the Relative Speed (v_rel) ---
    # The formula for relativistic velocity subtraction is: v_rel = (v1 - v2) / (1 - (v1 * v2) / c^2)
    # With c=1, this simplifies to: v_rel = (v1 - v2) / (1 - v1 * v2)
    try:
        calculated_v_rel = (v1 - v2) / (1 - v1 * v2)
    except ZeroDivisionError:
        return "Error: Division by zero in relative velocity calculation."

    # --- Step 2: Calculate the Total Energy (E_total) ---
    # The total energy is the sum of the individual relativistic energies: E_total = E1 + E2
    # The formula for relativistic energy is E = γ * m_rest * c^2, where γ = 1 / sqrt(1 - v^2/c^2)
    
    def lorentz_factor(v):
        if abs(v) >= 1:
            raise ValueError("Velocity cannot be greater than or equal to the speed of light.")
        return 1 / math.sqrt(1 - v**2)

    try:
        gamma1 = lorentz_factor(v1)
        gamma2 = lorentz_factor(v2)
    except ValueError as e:
        return f"Error in Lorentz factor calculation: {e}"

    # Calculate individual energies (in units of mc^2)
    E1 = gamma1 * m1_factor
    E2 = gamma2 * m2_factor
    
    # Calculate total energy
    calculated_E_total = E1 + E2

    # --- Step 3: Compare Calculated Values with the Answer ---
    # The values in the options are rounded to two decimal places.
    # We will compare our calculated values (also rounded) to the answer's values.
    
    # Check relative velocity
    # The calculated value is ~0.142857. Rounded to two decimal places, it is 0.14.
    if not math.isclose(round(calculated_v_rel, 2), answer_v_rel, rel_tol=1e-9):
        return (f"Incorrect. The relative speed calculation is wrong. "
                f"Calculated v_rel is approximately {calculated_v_rel:.4f}c, which rounds to {round(calculated_v_rel, 2)}c. "
                f"The answer C states v_rel = {answer_v_rel}c.")

    # Check total energy
    # The calculated value is ~5.9641. Rounded to two decimal places, it is 5.96.
    if not math.isclose(round(calculated_E_total, 2), answer_E_total, rel_tol=1e-9):
        return (f"Incorrect. The total energy calculation is wrong. "
                f"Calculated E_total is approximately {calculated_E_total:.4f}mc^2, which rounds to {round(calculated_E_total, 2)}mc^2. "
                f"The answer C states E_total = {answer_E_total}mc^2.")

    # If both checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_correctness_of_answer()
print(result)