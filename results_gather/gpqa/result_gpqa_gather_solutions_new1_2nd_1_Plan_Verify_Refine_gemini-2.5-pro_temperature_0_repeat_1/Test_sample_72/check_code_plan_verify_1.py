import math

def check_answer():
    """
    This function calculates the relative speed and total energy from first principles
    and compares them to the options provided in the question.
    """
    # Given values from the question
    m1_factor = 2.0
    v1_c = 0.6  # speed of astronaut 1 in units of c
    m2_factor = 3.0
    v2_c = 0.5  # speed of astronaut 2 in units of c

    # --- 1. Calculate Relative Speed (v_rel) ---
    # The formula for relativistic velocity subtraction for objects moving in the same direction is:
    # v_rel = (v1 - v2) / (1 - (v1 * v2) / c^2)
    # Since speeds are in units of c, the formula simplifies to:
    v_rel_calculated = (v1_c - v2_c) / (1 - v1_c * v2_c)

    # --- 2. Calculate Total Energy (E_total) ---
    # The formula for total relativistic energy is E = γmc², where γ = 1 / sqrt(1 - v²/c²)
    # The total energy of the system is E_total = E1 + E2

    # Calculate Lorentz factors
    try:
        gamma1 = 1 / math.sqrt(1 - v1_c**2)
        gamma2 = 1 / math.sqrt(1 - v2_c**2)
    except ValueError:
        return "Error: Calculation resulted in a value error, likely due to a speed >= c."

    # Calculate individual energies in units of mc^2
    e1_factor = gamma1 * m1_factor
    e2_factor = gamma2 * m2_factor

    # Calculate total energy in units of mc^2
    e_total_calculated = e1_factor + e2_factor

    # --- 3. Define the options from the question ---
    # A) v_rel =0.1c , E= 4.96 mc^2
    # B) v_rel = 0.14c , E= 5.96 mc^2
    # C) v_rel=0.14c, E=5mc^2
    # D) v_rel = 1.1c , E= mc^2
    
    # The correct option based on calculation should be B
    correct_option = {
        "v_rel": 0.14,
        "E_total": 5.96
    }

    # --- 4. Check correctness ---
    # We use math.isclose for floating point comparison, allowing for rounding differences.
    # The options are given to two decimal places, so a tolerance of 1e-2 is appropriate.
    v_rel_matches = math.isclose(v_rel_calculated, correct_option["v_rel"], rel_tol=1e-2)
    e_total_matches = math.isclose(e_total_calculated, correct_option["E_total"], rel_tol=1e-2)

    if v_rel_matches and e_total_matches:
        return "Correct"
    else:
        error_message = "Incorrect. The calculations do not match the proposed correct answer (Option B).\n"
        error_message += f"Calculated v_rel: {v_rel_calculated:.4f}c. Expected (rounded): {correct_option['v_rel']}c.\n"
        error_message += f"Calculated E_total: {e_total_calculated:.4f}mc^2. Expected (rounded): {correct_option['E_total']}mc^2.\n"
        if not v_rel_matches:
            error_message += "The relative velocity is incorrect.\n"
        if not e_total_matches:
            error_message += "The total energy is incorrect.\n"
        return error_message

# Run the check
result = check_answer()
print(result)