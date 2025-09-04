import math

def check_titration_answer():
    """
    This function calculates the pH at 25% titration and the equivalence point
    for the given problem and compares it to the provided answer 'A'.
    """
    # --- Problem Parameters ---
    acid_conc_initial = 0.05  # M
    acid_vol_initial_cm3 = 20.00 # cm3
    water_vol_cm3 = 20.00 # cm3
    base_conc = 0.1 # M
    Ka = 1.85e-5
    Kw = 1.0e-14

    # --- Proposed Answer (Option A) ---
    proposed_ph_25 = 4.26
    proposed_ph_eq = 8.52

    # --- Calculation for pH at 25% Titration ---
    # The Henderson-Hasselbalch equation is used: pH = pKa + log([A-]/[HA])
    # At 25% titration, 25% of the acid has been converted to its conjugate base.
    # The ratio of [A-]/[HA] is (0.25) / (1 - 0.25) = 0.25 / 0.75 = 1/3.
    # This ratio is independent of the initial concentration or dilution.
    pKa = -math.log10(Ka)
    calculated_ph_25 = pKa + math.log10(1/3)

    # --- Calculation for pH at the Equivalence Point ---
    # 1. Calculate initial moles of acetic acid
    acid_vol_initial_L = acid_vol_initial_cm3 / 1000.0
    initial_moles_acid = acid_conc_initial * acid_vol_initial_L

    # 2. Calculate the total volume of the acid solution after dilution
    total_acid_vol_L = (acid_vol_initial_cm3 + water_vol_cm3) / 1000.0

    # 3. Calculate the volume of NaOH needed to reach the equivalence point
    # At equivalence, moles of NaOH added = initial moles of acid
    vol_base_L = initial_moles_acid / base_conc

    # 4. Calculate the total volume of the solution at the equivalence point
    total_vol_at_eq_L = total_acid_vol_L + vol_base_L

    # 5. Calculate the concentration of the conjugate base (acetate, A-)
    # The moles of A- formed are equal to the initial moles of acid
    conc_A_minus_at_eq = initial_moles_acid / total_vol_at_eq_L

    # 6. The conjugate base hydrolyzes water: A- + H2O <=> HA + OH-
    # Calculate Kb for the conjugate base
    Kb = Kw / Ka

    # 7. Calculate the hydroxide ion concentration [OH-]
    # From the equilibrium expression Kb = [HA][OH-]/[A-], we can approximate
    # [OH-] = sqrt(Kb * [A-]), since Kb is small.
    conc_OH = math.sqrt(Kb * conc_A_minus_at_eq)

    # 8. Calculate pOH and then pH
    pOH = -math.log10(conc_OH)
    calculated_ph_eq = 14.0 - pOH

    # --- Verification Step ---
    # Compare the calculated values with the proposed answer, allowing for a small tolerance.
    tolerance = 0.01
    is_ph_25_correct = math.isclose(calculated_ph_25, proposed_ph_25, abs_tol=tolerance)
    is_ph_eq_correct = math.isclose(calculated_ph_eq, proposed_ph_eq, abs_tol=tolerance)

    if is_ph_25_correct and is_ph_eq_correct:
        return "Correct"
    else:
        error_messages = []
        if not is_ph_25_correct:
            error_messages.append(
                f"Constraint for pH at 25% titration is not satisfied. "
                f"Calculated pH is {calculated_ph_25:.2f}, but the answer provides {proposed_ph_25}."
            )
        if not is_ph_eq_correct:
            error_messages.append(
                f"Constraint for pH at the equivalence point is not satisfied. "
                f"Calculated pH is {calculated_ph_eq:.2f}, but the answer provides {proposed_ph_eq}."
            )
        return "\n".join(error_messages)

# Execute the check and print the result
result = check_titration_answer()
print(result)