import math

def check_titration_ph():
    """
    This function checks the correctness of the pH calculations for the given titration problem.
    It calculates the pH at 25% titration and at the equivalence point and compares them
    to the values in the proposed correct answer (Option C).
    """
    # --- Problem Parameters ---
    V_acid_initial_cm3 = 20.00
    M_acid_initial = 0.05  # Molarity of acetic acid
    V_water_cm3 = 20.00
    M_base = 0.1  # Molarity of NaOH
    Ka = 1.85e-5  # Ka of acetic acid
    Kw = 1.0e-14  # Ion product of water at 25°C

    # --- Step 1: Calculate the concentration of acetic acid after dilution ---
    # Convert volumes from cm³ to Liters
    V_acid_initial_L = V_acid_initial_cm3 / 1000
    V_water_L = V_water_cm3 / 1000

    # Calculate final volume and concentration after dilution
    V_acid_final_L = V_acid_initial_L + V_water_L
    # Using dilution formula M1V1 = M2V2 -> M2 = (M1V1)/V2
    M_acid_final = (M_acid_initial * V_acid_initial_L) / V_acid_final_L
    
    # Initial moles of acid in the diluted solution
    moles_acid = M_acid_final * V_acid_final_L

    # --- Step 2: Calculate pH at 25% titration ---
    # This is a buffer region, so we use the Henderson-Hasselbalch equation:
    # pH = pKa + log([A-]/[HA])
    pKa = -math.log10(Ka)
    
    # At 25% titration, the ratio of conjugate base [A-] to acid [HA] is 25/75 = 1/3
    ratio_A_HA = 25 / 75
    pH_at_25_percent = pKa + math.log10(ratio_A_HA)

    # --- Step 3: Calculate pH at the equivalence point ---
    # At the equivalence point, moles of base added = initial moles of acid.
    # All acid is converted to its conjugate base (acetate, A-).
    
    # Volume of NaOH needed to reach equivalence
    V_base_eq_L = moles_acid / M_base
    
    # Total volume of the solution at the equivalence point
    V_total_eq_L = V_acid_final_L + V_base_eq_L
    
    # Concentration of the conjugate base [A-] at the equivalence point
    conc_A_minus = moles_acid / V_total_eq_L
    
    # The pH is determined by the hydrolysis of the weak base A-:
    # A- + H2O <=> HA + OH-
    # We need the base dissociation constant, Kb
    Kb = Kw / Ka
    
    # From the equilibrium Kb = [HA][OH-]/[A-], we can approximate [OH-]
    # Kb ≈ [OH-]^2 / [A-], so [OH-] = sqrt(Kb * [A-])
    OH_concentration = math.sqrt(Kb * conc_A_minus)
    
    # Calculate pOH and then pH
    pOH = -math.log10(OH_concentration)
    pH_at_equivalence = 14.0 - pOH

    # --- Step 4: Check the correctness of the answer ---
    # The provided answer is C: 4.26; 8.52
    expected_pH_25 = 4.26
    expected_pH_eq = 8.52
    
    # We use a small tolerance for floating point comparisons
    tolerance = 0.01
    
    is_ph_25_correct = abs(pH_at_25_percent - expected_pH_25) < tolerance
    is_ph_eq_correct = abs(pH_at_equivalence - expected_pH_eq) < tolerance

    if is_ph_25_correct and is_ph_eq_correct:
        return "Correct"
    else:
        error_message = "Incorrect. The calculated values do not match the provided answer C (4.26; 8.52).\n"
        error_message += f"Calculated pH at 25% titration: {pH_at_25_percent:.2f}\n"
        error_message += f"Calculated pH at equivalence point: {pH_at_equivalence:.2f}\n"
        if not is_ph_25_correct:
            error_message += f"The pH at 25% titration is expected to be {expected_pH_25} but was calculated to be {pH_at_25_percent:.2f}.\n"
        if not is_ph_eq_correct:
            error_message += f"The pH at the equivalence point is expected to be {expected_pH_eq} but was calculated to be {pH_at_equivalence:.2f}.\n"
        return error_message

# Execute the check
result = check_titration_ph()
print(result)