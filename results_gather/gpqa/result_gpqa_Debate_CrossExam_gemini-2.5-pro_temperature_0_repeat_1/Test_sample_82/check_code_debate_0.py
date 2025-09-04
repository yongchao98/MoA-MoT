import math

def check_titration_ph_calculation():
    """
    This function checks the correctness of the given answer for a titration problem.
    It recalculates the pH at 25% titration and at the equivalence point based on the problem's parameters.
    """
    # --- Problem Parameters ---
    initial_conc_acid = 0.05  # M (mol/L)
    initial_vol_acid = 20.00  # cm^3
    vol_water_added = 20.00   # cm^3
    conc_NaOH = 0.1           # M (mol/L)
    Ka_acetic_acid = 1.85e-5
    Kw = 1.0e-14              # at 25 °C

    # --- Answer to be checked (Option C) ---
    expected_ph_25_percent = 4.26
    expected_ph_equiv = 8.52

    # --- Calculation for pH at 25% Titration ---

    # Step 1: Calculate the concentration of acetic acid after dilution.
    # M1V1 = M2V2 => M2 = (M1 * V1) / V2
    final_volume_after_dilution = initial_vol_acid + vol_water_added
    conc_acid_diluted = (initial_conc_acid * initial_vol_acid) / final_volume_after_dilution

    # Step 2: Use the Henderson-Hasselbalch equation for the buffer region.
    # pH = pKa + log([A-]/[HA])
    pKa = -math.log10(Ka_acetic_acid)
    
    # At 25% titration, 25% of the acid is converted to its conjugate base.
    # The ratio [A-]/[HA] is 25/75 or 1/3.
    ratio_A_minus_to_HA = 0.25 / 0.75
    calculated_ph_25_percent = pKa + math.log10(ratio_A_minus_to_HA)

    # --- Calculation for pH at the Equivalence Point ---

    # Step 1: Find the initial moles of diluted acid.
    # Moles = Molarity * Volume (in L)
    initial_moles_acid = conc_acid_diluted * (final_volume_after_dilution / 1000.0)

    # Step 2: Find the volume of NaOH needed to reach the equivalence point.
    # At equivalence, moles NaOH = moles acid.
    # Volume = Moles / Molarity
    vol_NaOH_equiv_L = initial_moles_acid / conc_NaOH
    vol_NaOH_equiv_cm3 = vol_NaOH_equiv_L * 1000.0

    # Step 3: Calculate the total volume at the equivalence point.
    total_volume_equiv_L = (final_volume_after_dilution + vol_NaOH_equiv_cm3) / 1000.0

    # Step 4: Calculate the concentration of the conjugate base (acetate) at the equivalence point.
    # Moles of acetate formed = initial moles of acid.
    conc_acetate_equiv = initial_moles_acid / total_volume_equiv_L

    # Step 5: Calculate the pH from the hydrolysis of the acetate ion.
    # A- + H2O <=> HA + OH-
    # Kb = Kw / Ka
    Kb = Kw / Ka_acetic_acid
    
    # [OH-] = sqrt(Kb * [A-])
    # This approximation is valid as Kb is small.
    conc_OH_minus = math.sqrt(Kb * conc_acetate_equiv)
    
    # pOH = -log10([OH-])
    pOH = -math.log10(conc_OH_minus)
    
    # pH = 14 - pOH
    calculated_ph_equiv = 14.0 - pOH

    # --- Verification ---
    # Check if the calculated values match the expected answer within a small tolerance.
    tolerance = 0.01

    is_ph_25_correct = abs(calculated_ph_25_percent - expected_ph_25_percent) < tolerance
    is_ph_equiv_correct = abs(calculated_ph_equiv - expected_ph_equiv) < tolerance

    if is_ph_25_correct and is_ph_equiv_correct:
        return "Correct"
    else:
        error_messages = []
        if not is_ph_25_correct:
            error_messages.append(
                f"Constraint for pH at 25% titration is not satisfied. "
                f"Expected value: {expected_ph_25_percent}, but calculated value is {calculated_ph_25_percent:.2f}."
            )
        if not is_ph_equiv_correct:
            error_messages.append(
                f"Constraint for pH at the equivalence point is not satisfied. "
                f"Expected value: {expected_ph_equiv}, but calculated value is {calculated_ph_equiv:.2f}."
            )
        return "\n".join(error_messages)

# Run the check
result = check_titration_ph_calculation()
print(result)