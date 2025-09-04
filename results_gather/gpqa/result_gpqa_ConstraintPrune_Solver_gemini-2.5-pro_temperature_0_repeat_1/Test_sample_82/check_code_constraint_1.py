import math

def check_titration_answer():
    """
    This function checks the correctness of the given answer for a titration problem.
    It recalculates the pH at 25% titration and at the equivalence point and
    compares the results with the provided answer.
    """
    # --- Problem Parameters ---
    # Initial acetic acid solution
    V_acid_initial_L = 20.00 / 1000  # 20.00 cm3 in Liters
    C_acid_initial = 0.05            # Molarity

    # Dilution
    V_water_L = 20.00 / 1000         # 20.00 cm3 in Liters

    # Titrant
    C_naoh = 0.1                     # Molarity of NaOH

    # Constants
    Ka_acid = 1.85e-5
    Kw = 1.0e-14
    pKa = -math.log10(Ka_acid)

    # --- Answer to be Checked ---
    # From option C
    expected_ph_25_percent = 4.26
    expected_ph_equivalence = 8.52

    # --- Calculation for pH at 25% Titration ---
    # At 25% titration, 25% of the acetic acid (HA) has been converted to its
    # conjugate base, acetate (A-), and 75% remains as HA.
    # This forms a buffer solution. We can use the Henderson-Hasselbalch equation:
    # pH = pKa + log([A-]/[HA])
    # The ratio of concentrations is the same as the ratio of moles (0.25 / 0.75).
    ratio_A_minus_to_HA = 0.25 / 0.75
    calculated_ph_25_percent = pKa + math.log10(ratio_A_minus_to_HA)

    # --- Calculation for pH at the Equivalence Point ---
    # 1. Calculate initial moles of acetic acid before dilution.
    moles_acid_initial = V_acid_initial_L * C_acid_initial

    # 2. At the equivalence point, moles of NaOH added equals initial moles of acid.
    #    Calculate the volume of NaOH required.
    V_naoh_at_equivalence_L = moles_acid_initial / C_naoh

    # 3. Calculate the total volume of the solution at the equivalence point.
    #    This is the volume of the diluted acid plus the volume of NaOH added.
    V_total_at_equivalence_L = (V_acid_initial_L + V_water_L) + V_naoh_at_equivalence_L

    # 4. At the equivalence point, all acetic acid has been converted to acetate (A-).
    #    The concentration of acetate is the initial moles of acid in the total new volume.
    C_acetate_at_equivalence = moles_acid_initial / V_total_at_equivalence_L

    # 5. Acetate hydrolyzes in water: A- + H2O <=> HA + OH-
    #    We need the base dissociation constant, Kb.
    Kb_acetate = Kw / Ka_acid

    # 6. Use the approximation to find [OH-]: Kb = [OH-]^2 / [A-]
    #    [OH-] = sqrt(Kb * [A-])
    concentration_OH = math.sqrt(Kb_acetate * C_acetate_at_equivalence)

    # 7. Calculate pOH and then pH.
    pOH = -math.log10(concentration_OH)
    calculated_ph_equivalence = 14.0 - pOH

    # --- Verification ---
    # Check if the calculated values are close to the expected values.
    # A tolerance is used to account for rounding in the problem's answer choices.
    tolerance = 0.015 
    
    is_ph_25_correct = math.isclose(calculated_ph_25_percent, expected_ph_25_percent, abs_tol=tolerance)
    is_ph_eq_correct = math.isclose(calculated_ph_equivalence, expected_ph_equivalence, abs_tol=tolerance)

    if is_ph_25_correct and is_ph_eq_correct:
        return "Correct"
    else:
        error_messages = []
        if not is_ph_25_correct:
            error_messages.append(
                f"pH at 25% titration is incorrect. Calculated value is {calculated_ph_25_percent:.2f}, but the answer is {expected_ph_25_percent}."
            )
        if not is_ph_eq_correct:
            error_messages.append(
                f"pH at the equivalence point is incorrect. Calculated value is {calculated_ph_equivalence:.2f}, but the answer is {expected_ph_equivalence}."
            )
        return "Incorrect. " + " ".join(error_messages)

# Run the check and print the result
result = check_titration_answer()
print(result)