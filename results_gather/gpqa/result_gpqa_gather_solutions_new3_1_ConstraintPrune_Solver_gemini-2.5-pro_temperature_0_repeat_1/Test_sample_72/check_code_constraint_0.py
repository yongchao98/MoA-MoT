import math

def check_correctness_of_answer():
    """
    This function calculates the physical quantities from the question and checks
    if they match the values provided in the selected answer option 'D'.
    """
    
    # --- Problem Parameters ---
    # Astronaut 1
    m1_factor = 2.0  # mass in units of m
    v1 = 0.6         # speed in units of c

    # Astronaut 2
    m2_factor = 3.0  # mass in units of m
    v2 = 0.5         # speed in units of c

    # --- Values from the selected answer (Option D) ---
    # D) v_rel = 0.14c , E= 5.96 mc^2
    expected_v_rel = 0.14  # in units of c
    expected_E_total = 5.96 # in units of mc^2

    # --- Calculation Step 1: Relative Speed (v_rel) ---
    # Using the relativistic velocity subtraction formula for co-linear motion:
    # v_rel = (v1 - v2) / (1 - (v1 * v2) / c^2)
    # Since speeds are in units of c, we can set c=1.
    try:
        calc_v_rel = (v1 - v2) / (1 - v1 * v2)
    except ZeroDivisionError:
        return "Error: Division by zero in relative speed calculation."

    # --- Calculation Step 2: Total Energy (E_total) ---
    # The total energy is the sum of individual relativistic energies:
    # E_total = E1 + E2 = (gamma1 * m1 * c^2) + (gamma2 * m2 * c^2)
    # We calculate the result in units of mc^2.
    
    try:
        # Lorentz factor (gamma) = 1 / sqrt(1 - v^2/c^2)
        gamma1 = 1 / math.sqrt(1 - v1**2)
        gamma2 = 1 / math.sqrt(1 - v2**2)

        # Individual energies in units of mc^2
        E1 = gamma1 * m1_factor
        E2 = gamma2 * m2_factor

        # Total energy
        calc_E_total = E1 + E2
    except ValueError:
        return "Error: Math domain error in energy calculation (e.g., speed >= c)."
    except ZeroDivisionError:
        return "Error: Division by zero in energy calculation (e.g., speed = c)."

    # --- Verification Step ---
    # Compare the calculated values with the expected values from option D.
    # We round the calculated values to two decimal places to match the precision of the options.
    
    # Check relative speed
    if round(calc_v_rel, 2) != expected_v_rel:
        return (f"Incorrect. The relative speed constraint is not satisfied. "
                f"Calculated relative speed is approximately {calc_v_rel:.4f}c, which rounds to {round(calc_v_rel, 2)}c. "
                f"The answer D states it is {expected_v_rel}c.")

    # Check total energy
    if round(calc_E_total, 2) != expected_E_total:
        return (f"Incorrect. The total energy constraint is not satisfied. "
                f"Calculated total energy is approximately {calc_E_total:.4f}mc^2, which rounds to {round(calc_E_total, 2)}mc^2. "
                f"The answer D states it is {expected_E_total}mc^2.")

    # If both values match the selected option D
    return "Correct"

# Execute the check and print the result
result = check_correctness_of_answer()
print(result)