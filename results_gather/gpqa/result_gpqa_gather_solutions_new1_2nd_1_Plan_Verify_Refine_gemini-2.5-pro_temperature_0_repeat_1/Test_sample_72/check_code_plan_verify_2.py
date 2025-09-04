import math

def check_correctness_of_answer():
    """
    This function verifies the final answer for the relativistic physics problem.
    It calculates the correct relative speed and total energy based on the principles
    of special relativity and compares them to the values in the selected answer option.
    """

    # --- Given parameters from the question ---
    # Astronaut 1
    m1_factor = 2.0  # Mass in units of 'm'
    v1 = 0.6       # Speed in units of 'c'

    # Astronaut 2
    m2_factor = 3.0  # Mass in units of 'm'
    v2 = 0.5       # Speed in units of 'c'

    # --- Step 1: Calculate the correct Relative Speed (v_rel) ---
    # The relativistic velocity subtraction formula for co-linear motion is used:
    # v_rel = (v1 - v2) / (1 - (v1 * v2) / c^2)
    # Since speeds are given as fractions of c, we can treat c as 1.
    try:
        calculated_v_rel = (v1 - v2) / (1 - v1 * v2)
    except ZeroDivisionError:
        return "Calculation Error: Division by zero occurred while calculating relative speed."

    # --- Step 2: Calculate the correct Total Energy (E_total) ---
    # The total energy of the system is the sum of the individual total relativistic energies.
    # E_total = E1 + E2 = (gamma1 * m1 * c^2) + (gamma2 * m2 * c^2)
    # The calculation will yield the total energy factor in units of mc^2.

    def lorentz_factor(v):
        """Calculates the Lorentz factor (gamma) for a given speed v (as a fraction of c)."""
        if abs(v) >= 1:
            raise ValueError("Speed cannot be greater than or equal to the speed of light.")
        return 1.0 / math.sqrt(1 - v**2)

    try:
        gamma1 = lorentz_factor(v1)
        gamma2 = lorentz_factor(v2)
    except ValueError as e:
        return f"Calculation Error: {e}"

    # Calculate the energy factor for each astronaut in units of mc^2
    E1_factor = gamma1 * m1_factor
    E2_factor = gamma2 * m2_factor

    # The total energy factor is the sum of the individual factors
    calculated_E_total_factor = E1_factor + E2_factor

    # --- Step 3: Parse the proposed answer from the prompt ---
    # The final answer given is <<<B>>>.
    # The options provided in the prompt are:
    # A) v_rel =0.1c , E= 4.96 mc^2
    # B) v_rel = 0.14c , E= 5.96 mc^2
    # C) v_rel=0.14c, E=5mc^2
    # D) v_rel = 1.1c , E= mc^2
    # We extract the numerical values from the chosen option B.
    proposed_v_rel = 0.14
    proposed_E_total_factor = 5.96

    # --- Step 4: Compare calculated values with the proposed answer's values ---
    # The options are given to two decimal places, so we round our calculated results
    # to the same precision for a direct comparison.
    
    rounded_calc_v_rel = round(calculated_v_rel, 2)
    rounded_calc_E_total_factor = round(calculated_E_total_factor, 2)

    # Check if both the relative speed and total energy match the proposed answer.
    v_rel_is_correct = (rounded_calc_v_rel == proposed_v_rel)
    E_total_is_correct = (rounded_calc_E_total_factor == proposed_E_total_factor)

    if v_rel_is_correct and E_total_is_correct:
        return "Correct"
    else:
        error_messages = []
        if not v_rel_is_correct:
            error_messages.append(
                f"The relative speed is incorrect. The calculated value is v_rel ≈ {calculated_v_rel:.4f}c, which rounds to {rounded_calc_v_rel}c, but the answer states it is {proposed_v_rel}c."
            )
        if not E_total_is_correct:
            error_messages.append(
                f"The total energy is incorrect. The calculated value is E ≈ {calculated_E_total_factor:.4f}mc^2, which rounds to {rounded_calc_E_total_factor}mc^2, but the answer states it is {proposed_E_total_factor}mc^2."
            )
        return "\n".join(error_messages)

# Execute the checking function and print the result.
result = check_correctness_of_answer()
print(result)