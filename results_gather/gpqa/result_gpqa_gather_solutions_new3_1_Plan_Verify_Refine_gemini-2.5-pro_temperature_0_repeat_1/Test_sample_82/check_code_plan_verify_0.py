import math

def check_titration_ph():
    """
    This function checks the correctness of the pH calculations for a weak acid-strong base titration.
    It recalculates the pH at 25% titration and at the equivalence point based on the problem's parameters.
    """
    # --- Problem Parameters ---
    initial_acid_vol_cm3 = 20.00
    initial_acid_conc_M = 0.05
    water_vol_cm3 = 20.00
    base_conc_M = 0.1
    Ka = 1.85e-5
    Kw = 1.0e-14
    temp_C = 25

    # --- Answer to be checked (Option B) ---
    # The provided answer is B, which corresponds to pH values of 4.26 and 8.52.
    expected_ph_25_percent = 4.26
    expected_ph_equiv = 8.52

    # --- Step 1: Calculate the concentration of acetic acid after dilution ---
    # M1V1 = M2V2 => M2 = (M1 * V1) / V2
    final_acid_vol_cm3 = initial_acid_vol_cm3 + water_vol_cm3
    diluted_acid_conc_M = (initial_acid_conc_M * initial_acid_vol_cm3) / final_acid_vol_cm3

    # --- Step 2: Calculate the pH at 25% titration ---
    # This is a buffer region, so we use the Henderson-Hasselbalch equation:
    # pH = pKa + log([A-]/[HA])
    pKa = -math.log10(Ka)
    # At 25% titration, 25% of the acid is converted to its conjugate base,
    # and 75% remains. The ratio [A-]/[HA] is 25/75 = 1/3.
    ratio_at_25_percent = 25.0 / 75.0
    calculated_ph_25_percent = pKa + math.log10(ratio_at_25_percent)

    # --- Step 3: Calculate the pH at the equivalence point ---
    # First, find the initial moles of acid in the diluted solution.
    final_acid_vol_L = final_acid_vol_cm3 / 1000.0
    initial_acid_moles = diluted_acid_conc_M * final_acid_vol_L

    # At equivalence, moles of base added = initial moles of acid.
    # Calculate the volume of base added.
    vol_base_added_L = initial_acid_moles / base_conc_M

    # Calculate the total volume at the equivalence point.
    total_vol_equiv_L = final_acid_vol_L + vol_base_added_L

    # At equivalence, all acid is converted to its conjugate base (A-).
    # The moles of A- are equal to the initial moles of acid.
    # Calculate the concentration of the conjugate base.
    conc_A_minus_equiv = initial_acid_moles / total_vol_equiv_L

    # The pH is determined by the hydrolysis of the conjugate base: A- + H2O <=> HA + OH-
    # Calculate Kb for the conjugate base.
    Kb = Kw / Ka

    # Use the approximation Kb = [OH-]^2 / [A-] to find [OH-].
    # This is valid since Kb is small.
    OH_conc = math.sqrt(Kb * conc_A_minus_equiv)

    # Calculate pOH and then pH.
    pOH = -math.log10(OH_conc)
    calculated_ph_equiv = 14.0 - pOH

    # --- Step 4: Verify the answer ---
    # Check if the calculated values match the expected answer within a small tolerance
    # to account for floating-point rounding.
    tolerance = 0.01
    
    is_ph_25_correct = abs(calculated_ph_25_percent - expected_ph_25_percent) < tolerance
    is_ph_equiv_correct = abs(calculated_ph_equiv - expected_ph_equiv) < tolerance

    if is_ph_25_correct and is_ph_equiv_correct:
        return "Correct"
    else:
        error_message = "Incorrect.\n"
        if not is_ph_25_correct:
            error_message += f"The pH at 25% titration is incorrect. Expected: {expected_ph_25_percent}, but calculated: {calculated_ph_25_percent:.2f}.\n"
        if not is_ph_equiv_correct:
            error_message += f"The pH at the equivalence point is incorrect. Expected: {expected_ph_equiv}, but calculated: {calculated_ph_equiv:.2f}.\n"
        
        # Provide a breakdown of the correct calculation for clarity.
        error_message += "\n--- Correct Calculation Breakdown ---\n"
        error_message += f"1. Diluted Acetic Acid Concentration: {diluted_acid_conc_M:.4f} M\n"
        error_message += f"2. pKa of Acetic Acid: {pKa:.3f}\n"
        error_message += f"3. pH at 25% Titration (pKa + log(1/3)): {calculated_ph_25_percent:.2f}\n"
        error_message += f"4. Moles of Acetic Acid: {initial_acid_moles:.5f} mol\n"
        error_message += f"5. Volume of NaOH at Equivalence: {vol_base_added_L * 1000:.2f} cm³\n"
        error_message += f"6. Total Volume at Equivalence: {total_vol_equiv_L * 1000:.2f} cm³\n"
        error_message += f"7. Acetate Concentration at Equivalence: {conc_A_minus_equiv:.4f} M\n"
        error_message += f"8. Kb for Acetate: {Kb:.3e}\n"
        error_message += f"9. [OH-] at Equivalence: {OH_conc:.3e} M\n"
        error_message += f"10. pOH at Equivalence: {pOH:.2f}\n"
        error_message += f"11. pH at Equivalence (14 - pOH): {calculated_ph_equiv:.2f}\n"
        
        return error_message.strip()

# Run the check
result = check_titration_ph()
print(result)