import math

def check_titration_ph():
    """
    This function checks the pH values at 25% titration and the equivalence point
    for the given chemistry problem.
    """
    # --- Problem Constants ---
    initial_vol_acid_cm3 = 20.00
    initial_conc_acid = 0.05  # M
    vol_water_added_cm3 = 20.00
    conc_NaOH = 0.1  # M
    Ka_acetic_acid = 1.85e-5
    Kw = 1.0e-14

    # --- Expected Answer Values from Option A ---
    expected_ph_25_percent = 4.26
    expected_ph_equivalence = 8.52

    # --- Step 1: Calculate the concentration of acetic acid after dilution ---
    final_vol_acid_cm3 = initial_vol_acid_cm3 + vol_water_added_cm3
    # Using M1V1 = M2V2
    conc_acid_diluted = (initial_conc_acid * initial_vol_acid_cm3) / final_vol_acid_cm3
    
    # Convert volume to Liters for molar calculations
    vol_acid_diluted_L = final_vol_acid_cm3 / 1000.0

    # --- Step 2: Calculate pH at 25% titration ---
    # This is a buffer region, so we use the Henderson-Hasselbalch equation:
    # pH = pKa + log([A-]/[HA])
    pKa = -math.log10(Ka_acetic_acid)
    
    # At 25% titration, the ratio of [A-]/[HA] is 25/75
    ratio_base_acid = 25.0 / 75.0
    
    calculated_ph_25_percent = pKa + math.log10(ratio_base_acid)

    # --- Step 3: Calculate pH at the equivalence point ---
    # At the equivalence point, all acid has been converted to its conjugate base.
    
    # Moles of acid initially present
    moles_acid = conc_acid_diluted * vol_acid_diluted_L
    
    # Volume of NaOH needed to reach equivalence (moles_acid = moles_base)
    vol_NaOH_eq_L = moles_acid / conc_NaOH
    
    # Total volume at equivalence point
    total_vol_eq_L = vol_acid_diluted_L + vol_NaOH_eq_L
    
    # Concentration of the conjugate base (acetate) at equivalence
    # Moles of acetate formed = initial moles of acid
    conc_acetate = moles_acid / total_vol_eq_L
    
    # Acetate hydrolyzes: A- + H2O <=> HA + OH-
    # We need Kb for acetate
    Kb_acetate = Kw / Ka_acetic_acid
    
    # Solve for [OH-] using the approximation Kb = x^2 / [A-]
    # where x = [OH-]
    conc_OH_squared = Kb_acetate * conc_acetate
    conc_OH = math.sqrt(conc_OH_squared)
    
    # Calculate pOH and then pH
    pOH = -math.log10(conc_OH)
    calculated_ph_equivalence = 14.00 - pOH

    # --- Step 4: Compare calculated values with the expected answer ---
    # Use a tolerance for floating point comparison
    tolerance = 0.01

    ph_25_correct = abs(calculated_ph_25_percent - expected_ph_25_percent) < tolerance
    ph_eq_correct = abs(calculated_ph_equivalence - expected_ph_equivalence) < tolerance

    if ph_25_correct and ph_eq_correct:
        return "Correct"
    else:
        error_messages = []
        if not ph_25_correct:
            error_messages.append(
                f"The pH at 25% titration is incorrect. "
                f"Expected: {expected_ph_25_percent}, "
                f"Calculated: {calculated_ph_25_percent:.2f}."
            )
        if not ph_eq_correct:
            error_messages.append(
                f"The pH at the equivalence point is incorrect. "
                f"Expected: {expected_ph_equivalence}, "
                f"Calculated: {calculated_ph_equivalence:.2f}."
            )
        return "\n".join(error_messages)

# Run the check and print the result
result = check_titration_ph()
print(result)