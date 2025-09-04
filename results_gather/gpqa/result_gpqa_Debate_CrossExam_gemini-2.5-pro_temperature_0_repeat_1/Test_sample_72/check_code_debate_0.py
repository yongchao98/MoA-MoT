import math

def check_astronaut_physics():
    """
    This function calculates the relative speed and total energy for the given physics problem
    and compares them to the values in the proposed answer.
    """
    # --- Problem Parameters ---
    # Let c = 1 and m = 1 for simplicity in calculation.
    # The final results will be in units of c and mc^2.
    m1 = 2.0  # Mass of astronaut 1 in units of m
    v1 = 0.6  # Speed of astronaut 1 in units of c
    m2 = 3.0  # Mass of astronaut 2 in units of m
    v2 = 0.5  # Speed of astronaut 2 in units of c

    # --- Values from the Proposed Answer (Option A) ---
    v_rel_answer = 0.14
    E_total_answer = 5.96

    # --- 1. Calculation of Relative Speed (v_rel) ---
    # The relativistic velocity subtraction formula for objects moving in the same direction is:
    # v_rel = (v_a - v_b) / (1 - (v_a * v_b) / c^2)
    # With c=1, the formula simplifies to: v_rel = (v1 - v2) / (1 - v1*v2)
    try:
        calculated_v_rel = (v1 - v2) / (1 - v1 * v2)
    except ZeroDivisionError:
        return "Error: Division by zero in relative speed calculation."

    # --- 2. Calculation of Total Energy (E_total) ---
    # The total energy is the sum of the individual relativistic energies: E_total = E1 + E2.
    # The formula for relativistic energy of a particle is E = gamma * m * c^2,
    # where the Lorentz factor gamma = 1 / sqrt(1 - v^2/c^2).
    # With c=1, E = gamma * m and gamma = 1 / sqrt(1 - v^2).

    # Energy of astronaut 1
    try:
        gamma1 = 1 / math.sqrt(1 - v1**2)
        E1 = gamma1 * m1
    except (ValueError, ZeroDivisionError) as e:
        return f"Error calculating energy for astronaut 1: {e}"

    # Energy of astronaut 2
    try:
        gamma2 = 1 / math.sqrt(1 - v2**2)
        E2 = gamma2 * m2
    except (ValueError, ZeroDivisionError) as e:
        return f"Error calculating energy for astronaut 2: {e}"

    # Total energy
    calculated_E_total = E1 + E2

    # --- 3. Verification ---
    # The answer provides values rounded to two decimal places.
    # We check if our calculated values, when rounded, match the answer.
    # A small tolerance is used for floating-point comparison.
    tolerance = 0.005 # A tolerance of half the last significant digit.

    v_rel_correct = abs(calculated_v_rel - v_rel_answer) < tolerance
    E_total_correct = abs(calculated_E_total - E_total_answer) < tolerance

    if v_rel_correct and E_total_correct:
        return "Correct"
    else:
        error_messages = []
        if not v_rel_correct:
            error_messages.append(
                f"Constraint failed: Relative speed calculation is incorrect.\n"
                f"Expected: v_rel ≈ {v_rel_answer}c\n"
                f"Calculated: v_rel ≈ {calculated_v_rel:.4f}c"
            )
        if not E_total_correct:
            error_messages.append(
                f"Constraint failed: Total energy calculation is incorrect.\n"
                f"Expected: E ≈ {E_total_answer}mc^2\n"
                f"Calculated: E ≈ {calculated_E_total:.4f}mc^2"
            )
        return "\n".join(error_messages)

# Run the check and print the result
result = check_astronaut_physics()
print(result)