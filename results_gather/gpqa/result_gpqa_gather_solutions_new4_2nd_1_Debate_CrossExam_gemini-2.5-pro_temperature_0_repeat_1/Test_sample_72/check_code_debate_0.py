import math

def check_answer():
    """
    Checks the correctness of the given answer for the relativistic physics problem.
    """
    # Given values from the question
    m1_factor = 2.0
    m2_factor = 3.0
    v1_factor = 0.6  # as a fraction of c
    v2_factor = 0.5  # as a fraction of c

    # The final answer provided is 'C'.
    # Let's extract the values from option C as presented in the final analysis.
    # C) v_rel = 0.14c , E= 5.96 mc^2
    expected_v_rel = 0.14
    expected_E_total = 5.96

    # --- Step 1: Calculate the relative speed (v_rel) ---
    # Using the relativistic velocity subtraction formula: v_rel = (v1 - v2) / (1 - (v1*v2)/c^2)
    # Since v1 and v2 are given as fractions of c, the formula simplifies to:
    # v_rel_factor = (v1_factor - v2_factor) / (1 - v1_factor * v2_factor)
    try:
        calc_v_rel = (v1_factor - v2_factor) / (1 - v1_factor * v2_factor)
    except ZeroDivisionError:
        return "Error: Division by zero in relative speed calculation."

    # --- Step 2: Calculate the total energy (E_total) ---
    # The total energy is the sum of individual relativistic energies: E_total = E1 + E2
    # where E = gamma * m * c^2 and gamma = 1 / sqrt(1 - (v/c)^2)
    # We will calculate the coefficient of mc^2.

    # Astronaut 1
    try:
        gamma1 = 1 / math.sqrt(1 - v1_factor**2)
        E1_factor = gamma1 * m1_factor
    except (ValueError, ZeroDivisionError) as e:
        return f"Error calculating energy for astronaut 1: {e}"

    # Astronaut 2
    try:
        gamma2 = 1 / math.sqrt(1 - v2_factor**2)
        E2_factor = gamma2 * m2_factor
    except (ValueError, ZeroDivisionError) as e:
        return f"Error calculating energy for astronaut 2: {e}"

    # Total Energy
    calc_E_total = E1_factor + E2_factor

    # --- Step 3: Compare calculated values with the expected answer ---
    # We use a small tolerance for floating-point comparisons.
    # Let's check if the rounded values match the answer's precision.
    v_rel_matches = math.isclose(round(calc_v_rel, 2), expected_v_rel)
    E_total_matches = math.isclose(round(calc_E_total, 2), expected_E_total)

    if v_rel_matches and E_total_matches:
        return "Correct"
    else:
        error_messages = []
        if not v_rel_matches:
            error_messages.append(f"Relative speed constraint not satisfied. Calculated v_rel is approximately {calc_v_rel:.4f}c, which rounds to {round(calc_v_rel, 2)}c. The answer expects {expected_v_rel}c.")
        if not E_total_matches:
            error_messages.append(f"Total energy constraint not satisfied. Calculated E_total is approximately {calc_E_total:.4f}mc^2, which rounds to {round(calc_E_total, 2)}mc^2. The answer expects {expected_E_total}mc^2.")
        return "\n".join(error_messages)

# Execute the check
result = check_answer()
print(result)