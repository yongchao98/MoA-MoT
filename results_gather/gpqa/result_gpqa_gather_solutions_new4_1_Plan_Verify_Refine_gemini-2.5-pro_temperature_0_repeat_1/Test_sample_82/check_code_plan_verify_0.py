import math

def check_titration_ph():
    """
    This function checks the correctness of the pH calculations for the titration problem.
    It recalculates the pH at 25% titration and at the equivalence point and compares
    it to the values given in the selected answer option.
    """
    # --- Given constants and initial values ---
    initial_vol_acid_cm3 = 20.00
    initial_conc_acid = 0.05  # M
    vol_water_cm3 = 20.00
    conc_naoh = 0.1  # M
    Ka = 1.85e-5
    Kw = 1.0e-14
    
    # The answer to check is C, which corresponds to pH values of 4.26 and 8.52
    expected_ph_25_percent = 4.26
    expected_ph_equivalence = 8.52

    # --- Convert volumes to Liters for calculation ---
    initial_vol_acid_L = initial_vol_acid_cm3 / 1000.0
    vol_water_L = vol_water_cm3 / 1000.0

    # --- Step 1: Calculate the concentration of acetic acid after dilution ---
    final_vol_acid_solution_L = initial_vol_acid_L + vol_water_L
    # Using dilution formula M1V1 = M2V2
    diluted_conc_acid = (initial_conc_acid * initial_vol_acid_L) / final_vol_acid_solution_L
    
    # --- Step 2: Calculate pH at 25% titration (Buffer region) ---
    # Use the Henderson-Hasselbalch equation: pH = pKa + log([A-]/[HA])
    pKa = -math.log10(Ka)
    # At 25% titration, the ratio of [A-]/[HA] is 25/75 = 1/3
    ratio = 25.0 / 75.0
    calculated_ph_25_percent = pKa + math.log10(ratio)

    # --- Step 3: Calculate pH at the equivalence point (Hydrolysis) ---
    # Initial moles of the diluted acid
    moles_acid = diluted_conc_acid * final_vol_acid_solution_L
    
    # Volume of NaOH needed to reach equivalence point
    # Moles of NaOH added = moles of acid
    vol_naoh_eq_L = moles_acid / conc_naoh
    
    # Total volume at equivalence point
    total_vol_eq_L = final_vol_acid_solution_L + vol_naoh_eq_L
    
    # Concentration of the conjugate base (acetate) at equivalence point
    # Moles of acetate = initial moles of acid
    conc_acetate = moles_acid / total_vol_eq_L
    
    # Calculate Kb for acetate
    Kb = Kw / Ka
    
    # Calculate [OH-] from hydrolysis: Kb = [OH-]^2 / [A-]
    # We assume [OH-] is small, so [A-] at equilibrium is approx. initial [A-]
    conc_oh = math.sqrt(Kb * conc_acetate)
    
    # Calculate pOH and then pH
    pOH = -math.log10(conc_oh)
    calculated_ph_equivalence = 14.0 - pOH

    # --- Step 4: Check the calculated values against the expected answer ---
    # Round the calculated values to two decimal places for comparison
    calculated_ph_25_percent_rounded = round(calculated_ph_25_percent, 2)
    calculated_ph_equivalence_rounded = round(calculated_ph_equivalence, 2)

    if calculated_ph_25_percent_rounded != expected_ph_25_percent:
        return (f"Incorrect: The pH at 25% titration is incorrect. "
                f"Expected: {expected_ph_25_percent}, Calculated: {calculated_ph_25_percent_rounded}.")

    if calculated_ph_equivalence_rounded != expected_ph_equivalence:
        return (f"Incorrect: The pH at the equivalence point is incorrect. "
                f"Expected: {expected_ph_equivalence}, Calculated: {calculated_ph_equivalence_rounded}.")

    return "Correct"

# Run the check
result = check_titration_ph()
print(result)