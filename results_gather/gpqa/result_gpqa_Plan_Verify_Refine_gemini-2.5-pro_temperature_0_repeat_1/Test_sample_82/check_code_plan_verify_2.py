import math

def check_titration_answer():
    """
    This function checks the correctness of the pH calculations for the titration of acetic acid with NaOH.
    It calculates the pH at 25% titration and at the equivalence point and compares them to the provided answer.
    """
    # --- Problem Constants and Initial Conditions ---
    VOL_ACID_INITIAL_CM3 = 20.00
    CONC_ACID_INITIAL = 0.05  # M
    VOL_WATER_CM3 = 20.00
    CONC_NAOH = 0.1  # M
    KA = 1.85e-5
    KW = 1.0e-14
    
    # The answer provided implies option B: pH at 25% is 4.26, pH at equivalence is 8.52
    EXPECTED_PH_25 = 4.26
    EXPECTED_PH_EQ = 8.52
    TOLERANCE = 0.02  # Allow for small rounding differences in pH values

    # --- Step 1: Calculate conditions after dilution ---
    vol_acid_initial_L = VOL_ACID_INITIAL_CM3 / 1000.0
    vol_water_L = VOL_WATER_CM3 / 1000.0
    
    initial_moles_acid = vol_acid_initial_L * CONC_ACID_INITIAL
    
    total_vol_diluted_L = vol_acid_initial_L + vol_water_L
    conc_acid_diluted = initial_moles_acid / total_vol_diluted_L

    # --- Step 2: Calculate pH at 25% Titration ---
    # This is a buffer region. Use the Henderson-Hasselbalch equation.
    # pH = pKa + log([A-]/[HA])
    # At 25% titration, the ratio [A-]/[HA] is 0.25 / 0.75 = 1/3.
    pKa = -math.log10(KA)
    calculated_pH_25 = pKa + math.log10(0.25 / 0.75)

    # --- Step 3: Calculate pH at the Equivalence Point ---
    # At equivalence, moles of NaOH added = initial moles of acid.
    # All acid is converted to its conjugate base, acetate (A-).
    vol_naoh_eq_L = initial_moles_acid / CONC_NAOH
    
    # Calculate total volume and concentration of acetate at the equivalence point.
    total_vol_eq_L = total_vol_diluted_L + vol_naoh_eq_L
    conc_acetate_eq = initial_moles_acid / total_vol_eq_L
    
    # Acetate hydrolyzes: A- + H2O <=> HA + OH-
    # We need Kb for acetate.
    Kb = KW / KA
    
    # Use the approximation Kb = [OH-]^2 / [A-]
    # [OH-] = sqrt(Kb * [A-])
    hydroxide_conc = math.sqrt(Kb * conc_acetate_eq)
    
    pOH = -math.log10(hydroxide_conc)
    calculated_pH_eq = 14.0 - pOH

    # --- Step 4: Verify the answer ---
    is_ph_25_correct = abs(calculated_pH_25 - EXPECTED_PH_25) <= TOLERANCE
    is_ph_eq_correct = abs(calculated_pH_eq - EXPECTED_PH_EQ) <= TOLERANCE

    if is_ph_25_correct and is_ph_eq_correct:
        return "Correct"
    else:
        error_messages = []
        if not is_ph_25_correct:
            error_messages.append(f"Incorrect pH at 25% titration. Calculated: {calculated_pH_25:.2f}, Expected: {EXPECTED_PH_25}.")
        if not is_ph_eq_correct:
            error_messages.append(f"Incorrect pH at equivalence point. Calculated: {calculated_pH_eq:.2f}, Expected: {EXPECTED_PH_EQ}.")
        return " ".join(error_messages)

# Run the check and print the result
result = check_titration_answer()
print(result)