import math

def check_titration_answer():
    """
    This function checks the correctness of the pH calculations for a titration problem.
    It recalculates the pH at 25% titration and at the equivalence point based on the problem's parameters.
    """
    # --- Given Parameters ---
    V_acid_initial = 20.00 / 1000  # L
    M_acid_initial = 0.05  # M
    V_water = 20.00 / 1000  # L
    M_NaOH = 0.1  # M
    Ka_acid = 1.85e-5
    Kw = 1.0e-14

    # --- The answer to check ---
    # The provided answer is B, which corresponds to the values (4.26, 8.52)
    expected_ph_25_percent = 4.26
    expected_ph_equivalence = 8.52

    # --- Step 1: Calculate the concentration of acetic acid after dilution ---
    V_acid_diluted = V_acid_initial + V_water
    moles_acid = M_acid_initial * V_acid_initial
    M_acid_diluted = moles_acid / V_acid_diluted

    # --- Step 2: Calculate the pH at 25% titration (Buffer Region) ---
    # Using the Henderson-Hasselbalch equation: pH = pKa + log([A-]/[HA])
    pKa = -math.log10(Ka_acid)
    # At 25% titration, the ratio of [A-]/[HA] is 25/75 = 1/3
    ratio_base_acid = 1/3
    calculated_ph_25_percent = pKa + math.log10(ratio_base_acid)

    # --- Step 3: Calculate the pH at the equivalence point ---
    # At the equivalence point, moles of NaOH added equals initial moles of acid.
    V_NaOH_added = moles_acid / M_NaOH
    
    # Total volume at equivalence point
    V_total_equivalence = V_acid_diluted + V_NaOH_added
    
    # Concentration of the conjugate base (acetate, A-) at equivalence point
    # Moles of A- formed are equal to the initial moles of acid
    M_conjugate_base = moles_acid / V_total_equivalence
    
    # pH is determined by the hydrolysis of the conjugate base: A- + H2O <=> HA + OH-
    # We need Kb for the conjugate base
    Kb_conjugate_base = Kw / Ka_acid
    
    # Using the approximation Kb = [OH-]^2 / [A-], we can find [OH-]
    # [OH-] = sqrt(Kb * [A-])
    conc_OH = math.sqrt(Kb_conjugate_base * M_conjugate_base)
    
    # Calculate pOH and then pH
    pOH = -math.log10(conc_OH)
    calculated_ph_equivalence = 14.00 - pOH

    # --- Step 4: Compare calculated values with the expected answer ---
    # Use a tolerance for floating point comparison
    tolerance = 0.01
    
    is_ph_25_correct = abs(calculated_ph_25_percent - expected_ph_25_percent) < tolerance
    is_ph_eq_correct = abs(calculated_ph_equivalence - expected_ph_equivalence) < tolerance

    if is_ph_25_correct and is_ph_eq_correct:
        return "Correct"
    else:
        error_message = "Incorrect. The provided answer does not match the calculated values.\n"
        if not is_ph_25_correct:
            error_message += f"Constraint failure (pH at 25% titration): Calculated pH is {calculated_ph_25_percent:.2f}, but the answer provides {expected_ph_25_percent}.\n"
        if not is_ph_eq_correct:
            error_message += f"Constraint failure (pH at equivalence point): Calculated pH is {calculated_ph_equivalence:.2f}, but the answer provides {expected_ph_equivalence}.\n"
        return error_message.strip()

# Run the check and print the result
result = check_titration_answer()
print(result)