import math

def check_titration_answer():
    """
    This function calculates the pH at 25% titration and at the equivalence point
    for the given chemistry problem and checks if the provided answer is correct.
    """
    # Given constants and initial values
    Ka = 1.85e-5
    Kw = 1.0e-14
    M_acid_initial = 0.05  # M
    V_acid_initial = 20.00 / 1000  # L
    V_water = 20.00 / 1000  # L
    M_base = 0.1  # M

    # The answer to check is B, which corresponds to (4.26, 8.52)
    expected_ph_25 = 4.26
    expected_ph_eq = 8.52

    # --- Step 1: Calculate pH at 25% titration ---
    # This is a buffer region, use Henderson-Hasselbalch equation:
    # pH = pKa + log([A-]/[HA])
    pKa = -math.log10(Ka)
    
    # At 25% titration, the ratio of [A-]/[HA] is 25/75 = 1/3
    ratio = 25.0 / 75.0
    calculated_ph_25 = pKa + math.log10(ratio)

    # --- Step 2: Calculate pH at the equivalence point ---
    # First, calculate the concentration of acetic acid after dilution
    V_acid_diluted = V_acid_initial + V_water
    M_acid_diluted = (M_acid_initial * V_acid_initial) / V_acid_diluted
    
    # Moles of acid to be titrated
    moles_acid = M_acid_diluted * V_acid_diluted
    
    # Volume of base needed to reach equivalence point
    V_base_eq = moles_acid / M_base
    
    # Total volume at equivalence point
    V_total_eq = V_acid_diluted + V_base_eq
    
    # Concentration of the conjugate base (acetate) at equivalence point
    M_conj_base = moles_acid / V_total_eq
    
    # Calculate Kb for the conjugate base
    Kb = Kw / Ka
    
    # Calculate [OH-] from hydrolysis: Kb = [OH-]^2 / [A-]
    # Assuming [OH-] is small, so [A-] - [OH-] is approx [A-]
    conc_OH = math.sqrt(Kb * M_conj_base)
    
    # Calculate pOH and then pH
    pOH_eq = -math.log10(conc_OH)
    calculated_ph_eq = 14.00 - pOH_eq

    # --- Step 3: Check the correctness of the answer ---
    # Use a tolerance for floating point comparison (answers are to 2 decimal places)
    tolerance = 0.01

    is_ph_25_correct = math.isclose(calculated_ph_25, expected_ph_25, abs_tol=tolerance)
    is_ph_eq_correct = math.isclose(calculated_ph_eq, expected_ph_eq, abs_tol=tolerance)

    if is_ph_25_correct and is_ph_eq_correct:
        return "Correct"
    else:
        error_message = "Incorrect. "
        if not is_ph_25_correct:
            error_message += f"The calculated pH at 25% titration is {calculated_ph_25:.2f}, but the answer expects {expected_ph_25}. "
        if not is_ph_eq_correct:
            error_message += f"The calculated pH at the equivalence point is {calculated_ph_eq:.2f}, but the answer expects {expected_ph_eq}."
        return error_message.strip()

# Run the check
result = check_titration_answer()
print(result)