import math

def check_titration_ph():
    """
    This function checks the correctness of the pH calculations for a weak acid-strong base titration problem.
    It recalculates the pH at 25% titration and at the equivalence point based on the problem's parameters.
    """
    # --- Given parameters from the question ---
    initial_acid_conc = 0.05  # M
    initial_acid_vol_cm3 = 20.00  # cm3
    water_vol_cm3 = 20.00  # cm3
    base_conc = 0.1  # M
    Ka_acetic_acid = 1.85e-5
    Kw = 1.0e-14
    temperature = 25  # °C

    # --- Expected answer values from option A ---
    expected_ph_25_percent = 4.26
    expected_ph_equivalence = 8.52

    # --- Step 1: Calculate the concentration of acetic acid after dilution ---
    final_acid_vol_cm3 = initial_acid_vol_cm3 + water_vol_cm3
    # Using the dilution formula M1V1 = M2V2
    diluted_acid_conc = (initial_acid_conc * initial_acid_vol_cm3) / final_acid_vol_cm3
    
    # --- Step 2: Calculate the pH at 25% titration ---
    # This is a buffer region, so we use the Henderson-Hasselbalch equation:
    # pH = pKa + log([A-]/[HA])
    
    # Calculate pKa
    pKa = -math.log10(Ka_acetic_acid)
    
    # At 25% titration, the ratio of [A-]/[HA] is 25/75 or 1/3
    ratio_A_HA = 25.0 / 75.0
    
    # Calculate pH
    calculated_ph_25_percent = pKa + math.log10(ratio_A_HA)

    # --- Step 3: Calculate the pH at the equivalence point ---
    # At the equivalence point, all the weak acid has been converted to its conjugate base.
    
    # 3a. Find the initial moles of the diluted acid
    final_acid_vol_L = final_acid_vol_cm3 / 1000.0
    moles_acid = diluted_acid_conc * final_acid_vol_L
    
    # 3b. At equivalence, moles of base added = moles of acid. Calculate the volume of NaOH added.
    vol_base_added_L = moles_acid / base_conc
    vol_base_added_cm3 = vol_base_added_L * 1000.0
    
    # 3c. Calculate the total volume of the solution at the equivalence point.
    total_vol_eq_L = final_acid_vol_L + vol_base_added_L
    
    # 3d. The moles of acetate (A-) formed are equal to the initial moles of acid.
    # Calculate the concentration of the acetate ion.
    conc_acetate = moles_acid / total_vol_eq_L
    
    # 3e. The acetate ion hydrolyzes water. Calculate Kb for acetate.
    Kb = Kw / Ka_acetic_acid
    
    # 3f. Set up the equilibrium expression for hydrolysis: Kb = [OH-]^2 / [A-]
    # Solve for [OH-], assuming [OH-] is small compared to [A-].
    OH_conc = math.sqrt(Kb * conc_acetate)
    
    # 3g. Calculate pOH and then pH.
    pOH = -math.log10(OH_conc)
    calculated_ph_equivalence = 14.0 - pOH

    # --- Step 4: Compare the calculated values with the expected answer ---
    # Round the calculated values to two decimal places for comparison.
    rounded_ph_25 = round(calculated_ph_25_percent, 2)
    rounded_ph_eq = round(calculated_ph_equivalence, 2)

    # Check if the calculated values match the expected values.
    if rounded_ph_25 == expected_ph_25_percent and rounded_ph_eq == expected_ph_equivalence:
        return "Correct"
    else:
        error_message = []
        if rounded_ph_25 != expected_ph_25_percent:
            error_message.append(
                f"The pH at 25% titration is incorrect. "
                f"Expected: {expected_ph_25_percent}, Calculated: {rounded_ph_25}."
            )
        if rounded_ph_eq != expected_ph_equivalence:
            error_message.append(
                f"The pH at the equivalence point is incorrect. "
                f"Expected: {expected_ph_equivalence}, Calculated: {rounded_ph_eq}."
            )
        return "\n".join(error_message)

# Run the check
result = check_titration_ph()
print(result)