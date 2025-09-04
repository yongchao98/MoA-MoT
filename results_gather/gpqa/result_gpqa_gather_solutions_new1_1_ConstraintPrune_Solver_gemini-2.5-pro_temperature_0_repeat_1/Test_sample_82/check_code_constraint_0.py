import math

def check_titration_answer():
    """
    This function calculates the pH at 25% titration and at the equivalence point
    for the given chemistry problem and checks if the provided answer is correct.
    """
    # --- Problem Parameters ---
    initial_vol_acid_cm3 = 20.00
    initial_conc_acid_M = 0.05
    vol_water_cm3 = 20.00
    conc_NaOH_M = 0.1
    Ka_acetic_acid = 1.85e-5
    Kw = 1.0e-14
    temperature_C = 25

    # --- Step 1: Calculate the concentration of acetic acid after dilution ---
    # Convert volumes from cm^3 to Liters
    initial_vol_acid_L = initial_vol_acid_cm3 / 1000
    vol_water_L = vol_water_cm3 / 1000

    # Calculate initial moles of acetic acid
    initial_moles_acid = initial_conc_acid_M * initial_vol_acid_L

    # Calculate the new volume and concentration after adding water
    final_vol_diluted_acid_L = initial_vol_acid_L + vol_water_L
    conc_acid_after_dilution_M = initial_moles_acid / final_vol_diluted_acid_L

    # --- Step 2: Calculate the pH at 25% titration ---
    # At 25% titration, the solution is a buffer. We use the Henderson-Hasselbalch equation:
    # pH = pKa + log([A-]/[HA])
    
    # Calculate pKa
    pKa = -math.log10(Ka_acetic_acid)
    
    # At 25% titration, 25% of the acid has been converted to its conjugate base (A-),
    # and 75% remains as acid (HA). The ratio [A-]/[HA] is 25/75 = 1/3.
    ratio_base_acid = 1 / 3
    
    calculated_ph_25_percent = pKa + math.log10(ratio_base_acid)

    # --- Step 3: Calculate the pH at the equivalence point ---
    # At the equivalence point, all the acetic acid has been neutralized by NaOH.
    # Moles of NaOH added = initial moles of acid
    vol_NaOH_equiv_L = initial_moles_acid / conc_NaOH_M
    
    # The total volume of the solution is the sum of the diluted acid and the added NaOH
    total_vol_equiv_L = final_vol_diluted_acid_L + vol_NaOH_equiv_L
    
    # The concentration of the conjugate base (acetate, A-) is the initial moles of acid
    # divided by the total volume at the equivalence point.
    conc_acetate_equiv_M = initial_moles_acid / total_vol_equiv_L
    
    # The pH is determined by the hydrolysis of the acetate ion: A- + H2O <=> HA + OH-
    # We need the base dissociation constant, Kb, for acetate.
    Kb_acetate = Kw / Ka_acetic_acid
    
    # We set up an equilibrium expression to find [OH-]: Kb = [HA][OH-] / [A-]
    # Assuming [OH-] = x, this becomes Kb = x^2 / ([A-] - x).
    # Since Kb is very small, we can approximate [A-] - x ≈ [A-].
    # So, [OH-] = sqrt(Kb * [A-])
    conc_OH = math.sqrt(Kb_acetate * conc_acetate_equiv_M)
    
    # Calculate pOH and then pH
    pOH = -math.log10(conc_OH)
    calculated_ph_equiv = 14.0 - pOH

    # --- Step 4: Check the correctness of the provided answer ---
    # The final answer from the analysis is D, which corresponds to (4.26, 8.52).
    expected_ph_25_percent = 4.26
    expected_ph_equiv = 8.52

    # Round the calculated values to two decimal places for a fair comparison.
    rounded_ph_25 = round(calculated_ph_25_percent, 2)
    rounded_ph_equiv = round(calculated_ph_equiv, 2)

    # Check if the calculated values match the expected values.
    is_ph25_correct = (rounded_ph_25 == expected_ph_25_percent)
    is_ph_equiv_correct = (rounded_ph_equiv == expected_ph_equiv)

    if is_ph25_correct and is_ph_equiv_correct:
        return "Correct"
    else:
        error_messages = []
        if not is_ph25_correct:
            error_messages.append(f"The pH at 25% titration is incorrect. Expected: {expected_ph_25_percent}, but calculated: {rounded_ph_25}.")
        if not is_ph_equiv_correct:
            error_messages.append(f"The pH at the equivalence point is incorrect. Expected: {expected_ph_equiv}, but calculated: {rounded_ph_equiv}.")
        return "Incorrect. " + " ".join(error_messages)

# Run the check and print the result.
result = check_titration_answer()
print(result)