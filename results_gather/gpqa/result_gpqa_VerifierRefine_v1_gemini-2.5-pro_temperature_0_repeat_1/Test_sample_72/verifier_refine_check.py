import math

def check_astronaut_physics_answer():
    """
    This function checks the correctness of the given answer for the relativistic physics problem.
    It calculates the relative speed and total energy of the system and compares them
    to the proposed answer.
    """
    # Given values from the question
    m1_factor = 2.0  # for mass 2m
    m2_factor = 3.0  # for mass 3m
    v1_factor = 0.6  # for speed 0.6c
    v2_factor = 0.5  # for speed 0.5c
    
    # --- Part 1: Calculate Relative Speed (v_rel) ---
    # The formula for relative velocity between two objects moving in the same direction is:
    # v_rel = (v1 - v2) / (1 - (v1 * v2) / c^2)
    # We can work with factors of c, so c=1 in the calculation.
    v_rel_numerator = v1_factor - v2_factor
    v_rel_denominator = 1 - (v1_factor * v2_factor)
    calculated_v_rel = v_rel_numerator / v_rel_denominator
    
    # --- Part 2: Calculate Total Energy (E_total) ---
    # The total energy of a particle is E = gamma * m * c^2
    # where gamma = 1 / sqrt(1 - (v/c)^2)
    # The total energy of the system is E_total = E1 + E2
    
    # Astronaut 1
    gamma1 = 1 / math.sqrt(1 - v1_factor**2)
    # Energy E1 is gamma1 * m1 * c^2 = gamma1 * (2m) * c^2
    # The coefficient of mc^2 is gamma1 * 2
    energy1_factor = gamma1 * m1_factor
    
    # Astronaut 2
    gamma2 = 1 / math.sqrt(1 - v2_factor**2)
    # Energy E2 is gamma2 * m2 * c^2 = gamma2 * (3m) * c^2
    # The coefficient of mc^2 is gamma2 * 3
    energy2_factor = gamma2 * m2_factor
    
    # Total Energy
    calculated_energy_factor = energy1_factor + energy2_factor
    
    # --- Part 3: Compare with the given answer (Option B) ---
    # Option B provides v_rel = 0.14c and E = 5.96 mc^2
    expected_v_rel = 0.14
    expected_energy_factor = 5.96
    
    # Use math.isclose for robust floating-point comparison, allowing for rounding.
    # The provided answer is rounded to 2 decimal places, so we set a tolerance.
    v_rel_correct = math.isclose(calculated_v_rel, expected_v_rel, rel_tol=1e-2)
    energy_correct = math.isclose(calculated_energy_factor, expected_energy_factor, rel_tol=1e-3)
    
    error_messages = []
    if not v_rel_correct:
        error_messages.append(
            f"Relative speed is incorrect. Expected approximately {expected_v_rel}c, but calculated {calculated_v_rel:.4f}c (which is 1/7 c)."
        )
    if not energy_correct:
        error_messages.append(
            f"Total energy is incorrect. Expected approximately {expected_energy_factor}mc^2, but calculated {calculated_energy_factor:.4f}mc^2."
        )
        
    if not error_messages:
        return "Correct"
    else:
        return "\n".join(error_messages)

# Run the check
result = check_astronaut_physics_answer()
print(result)