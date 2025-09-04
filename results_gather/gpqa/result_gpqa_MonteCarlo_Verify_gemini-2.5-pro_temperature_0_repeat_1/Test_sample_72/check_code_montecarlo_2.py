import math

def check_relativity_problem():
    """
    Checks the correctness of the answer to the relativistic astronauts problem.
    """
    # Given values from the question
    # Let m=1 and c=1 for calculation, as the answer is in terms of m and c.
    m1 = 2.0  # mass of astronaut 1 is 2m
    m2 = 3.0  # mass of astronaut 2 is 3m
    v1 = 0.6  # speed of astronaut 1 is 0.6c
    v2 = 0.5  # speed of astronaut 2 is 0.5c
    c = 1.0

    # Values from the proposed answer (Option A)
    v_rel_ans = 0.14  # in units of c
    E_total_ans = 5.96 # in units of mc^2

    # --- 1. Calculate the relative speed ---
    # Using the relativistic velocity subtraction formula: v_rel = (v1 - v2) / (1 - v1*v2/c^2)
    # We take the absolute value for speed.
    v_rel_calc = abs((v1 - v2) / (1 - (v1 * v2) / c**2))

    # --- 2. Calculate the total energy ---
    # Energy of a particle is E = gamma * m * c^2
    # where gamma = 1 / sqrt(1 - v^2/c^2)
    # Total energy is E_total = E1 + E2

    # Calculate Lorentz factor (gamma) for astronaut 1
    gamma1 = 1 / math.sqrt(1 - v1**2 / c**2)
    # Calculate energy for astronaut 1 (in units of mc^2)
    E1 = gamma1 * m1

    # Calculate Lorentz factor (gamma) for astronaut 2
    gamma2 = 1 / math.sqrt(1 - v2**2 / c**2)
    # Calculate energy for astronaut 2 (in units of mc^2)
    E2 = gamma2 * m2

    # Calculate total energy of the system
    E_total_calc = E1 + E2

    # --- 3. Compare calculated values with the answer ---
    # Use a tolerance for floating-point comparisons, appropriate for the precision in the answer.
    # The answer is given to two decimal places, so a tolerance of 0.01 is reasonable.
    v_rel_correct = math.isclose(v_rel_calc, v_rel_ans, rel_tol=0.01, abs_tol=0.01)
    E_total_correct = math.isclose(E_total_calc, E_total_ans, rel_tol=0.01, abs_tol=0.01)

    if v_rel_correct and E_total_correct:
        return "Correct"
    else:
        error_message = []
        if not v_rel_correct:
            error_message.append(
                f"Relative speed is incorrect. "
                f"Expected: {v_rel_ans}c, Calculated: {v_rel_calc:.4f}c."
            )
        if not E_total_correct:
            error_message.append(
                f"Total energy is incorrect. "
                f"Expected: {E_total_ans}mc^2, Calculated: {E_total_calc:.4f}mc^2."
            )
        return "\n".join(error_message)

# Run the check and print the result
result = check_relativity_problem()
print(result)