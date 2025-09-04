import math

def check_titration_answer():
    """
    This function calculates the pH at 25% titration and at the equivalence point
    for the given chemistry problem and compares it to the expected answer.
    """
    # --- Problem Parameters from the question ---
    initial_vol_acid_cm3 = 20.00
    initial_conc_acid = 0.05  # M
    vol_water_cm3 = 20.00
    conc_naoh = 0.1  # M
    Ka = 1.85e-5
    Kw = 1.0e-14
    pKw = 14.0

    # --- Expected Answer to Verify ---
    # The final answer from the analysis is B) 4.26; 8.52
    expected_ph_25 = 4.26
    expected_ph_eq = 8.52

    # --- Step 1: Calculate the concentration of acetic acid after dilution ---
    initial_vol_acid_L = initial_vol_acid_cm3 / 1000.0
    initial_moles_acid = initial_conc_acid * initial_vol_acid_L

    final_vol_diluted_cm3 = initial_vol_acid_cm3 + vol_water_cm3
    final_vol_diluted_L = final_vol_diluted_cm3 / 1000.0
    conc_acid_diluted = initial_moles_acid / final_vol_diluted_L

    # --- Step 2: Calculate the pH at 25% titration (Buffer Region) ---
    # Use the Henderson-Hasselbalch equation: pH = pKa + log([A-]/[HA])
    pKa = -math.log10(Ka)
    # At 25% titration, 25% of the acid is converted to its conjugate base,
    # and 75% remains. The ratio [A-]/[HA] is 25/75 = 1/3.
    ratio_A_HA = 0.25 / 0.75
    ph_25_percent = pKa + math.log10(ratio_A_HA)

    # --- Step 3: Calculate the pH at the equivalence point (Hydrolysis) ---
    # At the equivalence point, moles of NaOH added = initial moles of acid.
    # Calculate the volume of NaOH needed to reach equivalence.
    vol_naoh_eq_L = initial_moles_acid / conc_naoh
    
    # Calculate the total volume at the equivalence point.
    total_vol_eq_L = final_vol_diluted_L + vol_naoh_eq_L
    
    # Calculate the concentration of the conjugate base (acetate) at this point.
    conc_acetate_eq = initial_moles_acid / total_vol_eq_L
    
    # The acetate ion hydrolyzes: A- + H2O <=> HA + OH-
    # Calculate Kb for the acetate ion.
    Kb = Kw / Ka
    
    # Solve for [OH-] using the approximation Kb = [OH-]^2 / [A-],
    # which is valid because Kb is very small.
    # x^2 = Kb * [A-]
    conc_oh_squared = Kb * conc_acetate_eq
    conc_oh = math.sqrt(conc_oh_squared)
    
    # Calculate pOH and then convert to pH.
    pOH_eq = -math.log10(conc_oh)
    pH_eq = pKw - pOH_eq

    # --- Step 4: Verification ---
    # Round the calculated values to two decimal places for comparison.
    calculated_ph_25_rounded = round(ph_25_percent, 2)
    calculated_ph_eq_rounded = round(pH_eq, 2)

    # Check if the calculated values match the expected answer within a small tolerance.
    tolerance = 0.01
    is_ph25_correct = abs(calculated_ph_25_rounded - expected_ph_25) <= tolerance
    is_ph_eq_correct = abs(calculated_ph_eq_rounded - expected_ph_eq) <= tolerance

    if is_ph25_correct and is_ph_eq_correct:
        return "Correct"
    else:
        reasons = []
        if not is_ph25_correct:
            reasons.append(
                f"The pH at 25% titration is incorrect. "
                f"Expected: {expected_ph_25}, Calculated: {calculated_ph_25_rounded:.2f}."
            )
        if not is_ph_eq_correct:
            reasons.append(
                f"The pH at the equivalence point is incorrect. "
                f"Expected: {expected_ph_eq}, Calculated: {calculated_ph_eq_rounded:.2f}."
            )
        
        # Provide detailed calculation steps in the error message for clarity.
        details = (
            f"\n--- Calculation Details ---\n"
            f"Initial moles of acid: {initial_moles_acid:.5f} mol\n"
            f"Diluted acid concentration: {conc_acid_diluted:.4f} M\n"
            f"pKa: {-math.log10(Ka):.3f}\n"
            f"Calculated pH at 25% titration: {ph_25_percent:.3f}\n"
            f"Volume of NaOH at equivalence: {vol_naoh_eq_L*1000:.2f} cm3\n"
            f"Concentration of acetate at equivalence: {conc_acetate_eq:.4f} M\n"
            f"Kb for acetate: {Kb:.3e}\n"
            f"[OH-] at equivalence: {conc_oh:.3e} M\n"
            f"pOH at equivalence: {pOH_eq:.3f}\n"
            f"Calculated pH at equivalence: {pH_eq:.3f}\n"
        )
        
        return "Incorrect. " + " ".join(reasons) + details

# Run the check and print the result.
result = check_titration_answer()
print(result)