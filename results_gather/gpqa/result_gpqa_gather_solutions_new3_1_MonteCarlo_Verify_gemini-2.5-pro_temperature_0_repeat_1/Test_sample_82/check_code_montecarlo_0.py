import math

def check_titration_answer():
    """
    This function calculates the pH at 25% titration and at the equivalence point
    for the given chemistry problem and checks if it matches the provided answer.
    """
    # --- Problem Parameters ---
    initial_acid_vol_cm3 = 20.00
    initial_acid_conc_M = 0.05
    water_vol_cm3 = 20.00
    base_conc_M = 0.1
    Ka = 1.85e-5
    Kw = 1.0e-14
    
    # The answer from the LLM to be checked (Option B)
    # pH at 25% titration; pH at equivalence point
    expected_ph_values = (4.26, 8.52)

    # --- Step 1: Dilution Calculation ---
    # Calculate the initial moles of acetic acid before dilution
    initial_moles_acid = (initial_acid_vol_cm3 / 1000.0) * initial_acid_conc_M
    
    # Calculate the volume of the acid solution after dilution
    diluted_acid_vol_cm3 = initial_acid_vol_cm3 + water_vol_cm3
    diluted_acid_vol_L = diluted_acid_vol_cm3 / 1000.0
    
    # Calculate the concentration of the diluted acetic acid
    diluted_acid_conc_M = initial_moles_acid / diluted_acid_vol_L

    # --- Step 2: pH at 25% Titration ---
    # At 25% titration, the solution is a buffer.
    # We use the Henderson-Hasselbalch equation: pH = pKa + log([A-]/[HA])
    pKa = -math.log10(Ka)
    
    # At 25% titration, the ratio of [A-]/[HA] is 25/75 or 1/3.
    ratio_A_HA = 25.0 / 75.0
    calculated_ph_25_percent = pKa + math.log10(ratio_A_HA)

    # --- Step 3: pH at the Equivalence Point ---
    # At the equivalence point, moles of base added = initial moles of acid.
    # Calculate the volume of NaOH added to reach the equivalence point.
    vol_base_added_L = initial_moles_acid / base_conc_M
    
    # Calculate the total volume of the solution at the equivalence point.
    total_vol_equiv_L = diluted_acid_vol_L + vol_base_added_L
    
    # At the equivalence point, all acid is converted to its conjugate base (acetate, A-).
    # The moles of acetate are equal to the initial moles of acid.
    # Calculate the concentration of the acetate ion.
    conc_A_minus_equiv = initial_moles_acid / total_vol_equiv_L
    
    # The acetate ion hydrolyzes, making the solution basic.
    # A- + H2O <=> HA + OH-
    # We need the base dissociation constant, Kb.
    Kb = Kw / Ka
    
    # Calculate the hydroxide ion concentration [OH-] using the approximation:
    # Kb = [HA][OH-]/[A-] ≈ [OH-]^2 / [A-]
    conc_OH = math.sqrt(Kb * conc_A_minus_equiv)
    
    # Calculate pOH and then pH.
    pOH = -math.log10(conc_OH)
    calculated_ph_equiv = 14.0 - pOH

    # --- Step 4: Verification ---
    # Compare the calculated values with the expected values from the answer.
    # We use a tolerance to account for floating-point arithmetic and rounding.
    is_ph_25_correct = math.isclose(calculated_ph_25_percent, expected_ph_values[0], abs_tol=0.01)
    is_ph_equiv_correct = math.isclose(calculated_ph_equiv, expected_ph_values[1], abs_tol=0.01)

    if is_ph_25_correct and is_ph_equiv_correct:
        return "Correct"
    else:
        error_messages = []
        if not is_ph_25_correct:
            error_messages.append(
                f"Constraint for pH at 25% titration is not satisfied. "
                f"Calculated pH is {calculated_ph_25_percent:.2f}, but the answer provides {expected_ph_values[0]}."
            )
        if not is_ph_equiv_correct:
            error_messages.append(
                f"Constraint for pH at the equivalence point is not satisfied. "
                f"Calculated pH is {calculated_ph_equiv:.2f}, but the answer provides {expected_ph_values[1]}."
            )
        return "\n".join(error_messages)

# Execute the check and print the result
result = check_titration_answer()
print(result)