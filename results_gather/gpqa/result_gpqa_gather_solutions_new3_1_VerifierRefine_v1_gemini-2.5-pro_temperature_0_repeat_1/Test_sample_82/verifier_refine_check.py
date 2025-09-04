import math

def check_titration_ph():
    """
    This function calculates the pH at 25% titration and at the equivalence point
    for the given chemistry problem and checks if the provided answer is correct.
    """
    # --- Given constants and initial values ---
    initial_acid_vol_cm3 = 20.00
    initial_acid_conc = 0.05  # M
    water_vol_cm3 = 20.00
    naoh_conc = 0.1  # M
    Ka = 1.85e-5
    Kw = 1.0e-14
    pKw = 14.0

    # --- Step 1: Calculate the concentration of acetic acid after dilution ---
    # Convert volumes from cm³ to L
    initial_acid_vol_L = initial_acid_vol_cm3 / 1000.0
    water_vol_L = water_vol_cm3 / 1000.0

    # Calculate initial moles of acetic acid
    initial_moles_acid = initial_acid_conc * initial_acid_vol_L

    # Calculate the final volume after dilution
    diluted_acid_vol_L = initial_acid_vol_L + water_vol_L

    # Calculate the concentration of the diluted acetic acid
    diluted_acid_conc = initial_moles_acid / diluted_acid_vol_L

    # --- Step 2: Calculate the pH at 25% titration ---
    # This is a buffer region, so we use the Henderson-Hasselbalch equation:
    # pH = pKa + log([A-]/[HA])
    
    # Calculate pKa
    pKa = -math.log10(Ka)

    # At 25% titration, the ratio of [A-]/[HA] is 25/75 or 1/3
    ratio_A_HA = 25.0 / 75.0
    
    # Calculate pH
    ph_at_25_percent = pKa + math.log10(ratio_A_HA)

    # --- Step 3: Calculate the pH at the equivalence point ---
    # At the equivalence point, all acid has been converted to its conjugate base (acetate, A-).
    
    # Moles of A- formed are equal to the initial moles of acid
    moles_A_minus = initial_moles_acid

    # Calculate the volume of NaOH needed to reach the equivalence point
    vol_naoh_equiv_L = initial_moles_acid / naoh_conc

    # Calculate the total volume at the equivalence point
    total_vol_equiv_L = diluted_acid_vol_L + vol_naoh_equiv_L

    # Calculate the concentration of the conjugate base [A-] at the equivalence point
    conc_A_minus = moles_A_minus / total_vol_equiv_L

    # The pH is determined by the hydrolysis of A-: A- + H2O <=> HA + OH-
    # We need the base dissociation constant, Kb
    Kb = Kw / Ka

    # Set up the equilibrium expression: Kb = [HA][OH-] / [A-]
    # Let x = [OH-]. Then Kb = x^2 / ([A-] - x).
    # Since Kb is very small, we can approximate [A-] - x ≈ [A-].
    # So, x^2 = Kb * [A-]
    conc_OH_squared = Kb * conc_A_minus
    conc_OH = math.sqrt(conc_OH_squared)

    # Calculate pOH and then pH
    pOH = -math.log10(conc_OH)
    ph_at_equivalence = pKw - pOH

    # --- Step 4: Check the correctness of the provided answer ---
    # The provided answer is D, which corresponds to pH values of 4.26 and 8.52.
    expected_ph_25 = 4.26
    expected_ph_equiv = 8.52
    
    # Check if the calculated values match the expected values within a small tolerance
    tolerance = 0.01
    is_ph_25_correct = math.isclose(ph_at_25_percent, expected_ph_25, rel_tol=tolerance)
    is_ph_equiv_correct = math.isclose(ph_at_equivalence, expected_ph_equiv, rel_tol=tolerance)

    if is_ph_25_correct and is_ph_equiv_correct:
        return "Correct"
    else:
        reasons = []
        if not is_ph_25_correct:
            reasons.append(f"Incorrect pH at 25% titration. Expected: {expected_ph_25}, Calculated: {ph_at_25_percent:.2f}")
        if not is_ph_equiv_correct:
            reasons.append(f"Incorrect pH at equivalence point. Expected: {expected_ph_equiv}, Calculated: {ph_at_equivalence:.2f}")
        return "\n".join(reasons)

# Run the check
result = check_titration_ph()
print(result)