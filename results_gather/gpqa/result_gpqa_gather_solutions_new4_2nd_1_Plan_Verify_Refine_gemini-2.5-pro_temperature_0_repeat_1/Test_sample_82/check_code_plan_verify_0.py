import math

def check_titration_ph():
    """
    This function calculates the pH at 25% titration and at the equivalence point
    for the given problem and checks if the provided answer is correct.
    """
    # --- Problem Parameters ---
    initial_conc_acid = 0.05  # M (mol/L)
    initial_vol_acid = 20.00 / 1000  # L (from cm3)
    vol_water = 20.00 / 1000         # L (from cm3)
    conc_naoh = 0.1                  # M (mol/L)
    Ka_acid = 1.85e-5
    Kw = 1.0e-14
    
    # The provided answer corresponds to option C
    expected_ph_25_percent = 4.26
    expected_ph_equivalence = 8.52

    # --- Step 1: Calculate concentration of acetic acid after dilution ---
    # Using the dilution formula M1V1 = M2V2
    diluted_vol_acid = initial_vol_acid + vol_water
    try:
        diluted_conc_acid = (initial_conc_acid * initial_vol_acid) / diluted_vol_acid
    except ZeroDivisionError:
        return "Error: Total volume of diluted acid is zero."

    # --- Step 2: Calculate pH at 25% titration ---
    # This is a buffer region, so we use the Henderson-Hasselbalch equation:
    # pH = pKa + log([A-]/[HA])
    try:
        pKa = -math.log10(Ka_acid)
    except ValueError:
        return "Error: Ka must be positive to calculate pKa."
        
    # At 25% titration, the ratio of [A-]/[HA] is 25/75 = 1/3
    ratio_base_acid = 1/3
    calculated_ph_25_percent = pKa + math.log10(ratio_base_acid)

    # --- Step 3: Calculate pH at the equivalence point ---
    # At the equivalence point, all acid is converted to its conjugate base (acetate).
    # First, find the moles of acid being titrated.
    moles_acid = diluted_conc_acid * diluted_vol_acid
    
    # Next, find the volume of NaOH needed to reach equivalence.
    # Moles of NaOH = Moles of acid
    try:
        vol_naoh_at_equivalence = moles_acid / conc_naoh
    except ZeroDivisionError:
        return "Error: Concentration of NaOH cannot be zero."
        
    # Calculate the total volume at the equivalence point.
    total_vol_at_equivalence = diluted_vol_acid + vol_naoh_at_equivalence
    
    # Calculate the concentration of the conjugate base (acetate) at this point.
    try:
        conc_acetate_at_equivalence = moles_acid / total_vol_at_equivalence
    except ZeroDivisionError:
        return "Error: Total volume at equivalence point is zero."

    # Now, calculate the pH from the hydrolysis of the acetate ion.
    # A- + H2O <=> HA + OH-
    # We need Kb for the acetate ion.
    try:
        Kb_acetate = Kw / Ka_acid
    except ZeroDivisionError:
        return "Error: Ka cannot be zero."
        
    # Use the approximation Kb = [OH-]^2 / [A-] to find [OH-].
    # This assumes that the concentration of OH- from hydrolysis is small.
    conc_OH = math.sqrt(Kb_acetate * conc_acetate_at_equivalence)
    
    # Calculate pOH and then pH.
    pOH_at_equivalence = -math.log10(conc_OH)
    calculated_ph_equivalence = 14.0 - pOH_at_equivalence

    # --- Step 4: Check the correctness of the answer ---
    # Compare calculated values with the expected values, allowing for a small tolerance for rounding.
    tolerance = 0.01
    
    ph_25_correct = math.isclose(calculated_ph_25_percent, expected_ph_25_percent, rel_tol=tolerance)
    ph_equiv_correct = math.isclose(calculated_ph_equivalence, expected_ph_equivalence, rel_tol=tolerance)

    if ph_25_correct and ph_equiv_correct:
        return "Correct"
    else:
        error_message = "Incorrect. The calculated values do not match the provided answer.\n"
        if not ph_25_correct:
            error_message += f"Constraint failed for pH at 25% titration.\n"
            error_message += f"Expected: {expected_ph_25_percent}, Calculated: {calculated_ph_25_percent:.2f}\n"
        if not ph_equiv_correct:
            error_message += f"Constraint failed for pH at the equivalence point.\n"
            error_message += f"Expected: {expected_ph_equivalence}, Calculated: {calculated_ph_equivalence:.2f}\n"
        return error_message

# Execute the check
result = check_titration_ph()
print(result)