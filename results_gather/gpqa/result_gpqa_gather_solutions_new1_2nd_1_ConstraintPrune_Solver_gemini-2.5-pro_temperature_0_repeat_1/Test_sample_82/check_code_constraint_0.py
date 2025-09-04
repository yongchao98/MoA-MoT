import math

def check_answer():
    """
    This function checks the correctness of the given answer for the titration problem.
    It recalculates the pH values at 25% titration and the equivalence point based on the problem statement.
    """
    # --- Given parameters from the question ---
    initial_acid_vol_cm3 = 20.00
    initial_acid_conc = 0.05  # M
    water_vol_cm3 = 20.00
    base_conc = 0.1  # M
    Ka_acid = 1.85e-5
    Kw = 1.0e-14

    # --- The answer to check (Option A) ---
    # The final answer provided is <<<A>>>, which corresponds to the pair (4.26, 8.52)
    expected_ph_25_percent = 4.26
    expected_ph_equivalence = 8.52

    # --- Step 1: Calculate the concentration of acetic acid after dilution ---
    final_acid_vol_cm3 = initial_acid_vol_cm3 + water_vol_cm3
    # Using the dilution formula M1V1 = M2V2
    diluted_acid_conc = (initial_acid_conc * initial_acid_vol_cm3) / final_acid_vol_cm3

    # --- Step 2: Calculate the pH at 25% titration ---
    # At this point, the solution is a buffer. Use the Henderson-Hasselbalch equation:
    # pH = pKa + log([A-]/[HA])
    pKa = -math.log10(Ka_acid)
    # At 25% titration, the ratio of [conjugate base]/[acid] is 25/75 = 1/3
    ratio_base_acid = 25 / 75
    calculated_ph_25_percent = pKa + math.log10(ratio_base_acid)

    # --- Step 3: Calculate the pH at the equivalence point ---
    # First, find the initial moles of the diluted acid
    initial_moles_acid = diluted_acid_conc * (final_acid_vol_cm3 / 1000)

    # At the equivalence point, moles of base added equals initial moles of acid.
    # Calculate the volume of base added.
    vol_base_added_L = initial_moles_acid / base_conc

    # Calculate the total volume of the solution at the equivalence point.
    total_vol_L = (final_acid_vol_cm3 / 1000) + vol_base_added_L

    # At the equivalence point, all acid is converted to its conjugate base.
    # The moles of conjugate base equal the initial moles of acid.
    # Calculate the concentration of the conjugate base.
    conc_conjugate_base = initial_moles_acid / total_vol_L

    # The pH is determined by the hydrolysis of the conjugate base: A- + H2O <=> HA + OH-
    # We need the base dissociation constant, Kb.
    Kb = Kw / Ka_acid

    # Use the equilibrium expression for hydrolysis: Kb = [HA][OH-]/[A-]
    # We can approximate [OH-] = sqrt(Kb * [A-]) since Kb is small.
    conc_OH = math.sqrt(Kb * conc_conjugate_base)

    # Calculate pOH and then pH.
    pOH = -math.log10(conc_OH)
    calculated_ph_equivalence = 14 - pOH

    # --- Step 4: Compare calculated values with the expected answer ---
    # Round the calculated values to two decimal places for a fair comparison.
    rounded_ph_25 = round(calculated_ph_25_percent, 2)
    rounded_ph_eq = round(calculated_ph_equivalence, 2)

    # Check for correctness
    is_ph_25_correct = math.isclose(rounded_ph_25, expected_ph_25_percent, abs_tol=0.01)
    is_ph_eq_correct = math.isclose(rounded_ph_eq, expected_ph_equivalence, abs_tol=0.01)

    if is_ph_25_correct and is_ph_eq_correct:
        return "Correct"
    else:
        error_messages = []
        if not is_ph_25_correct:
            error_messages.append(
                f"The pH at 25% titration is incorrect. The provided answer is {expected_ph_25_percent}, but the calculated value is {rounded_ph_25}."
            )
        if not is_ph_eq_correct:
            error_messages.append(
                f"The pH at the equivalence point is incorrect. The provided answer is {expected_ph_equivalence}, but the calculated value is {rounded_ph_eq}."
            )
        return "\n".join(error_messages)

# Execute the check and print the result
result = check_answer()
print(result)