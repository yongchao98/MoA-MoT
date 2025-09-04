import math

def check_answer():
    """
    This function calculates the pH at 25% titration and at the equivalence point
    for the given problem and compares it to the provided answer (Option C).
    """
    # --- Problem Constraints & Given Values ---
    V_acid_initial_cm3 = 20.00
    C_acid_initial_M = 0.05
    V_water_cm3 = 20.00
    C_base_M = 0.1
    Ka_acid = 1.85e-5
    Kw = 1.0e-14

    # The answer to be checked (Option C)
    expected_ph_at_25_percent = 4.26
    expected_ph_at_equivalence = 8.52

    # --- Calculation for pH at 25% Titration ---

    # 1. Calculate pKa of acetic acid
    pKa = -math.log10(Ka_acid)

    # 2. Use the Henderson-Hasselbalch equation: pH = pKa + log([A-]/[HA])
    # At 25% titration, 25% of the acid (HA) has been converted to its
    # conjugate base (A-). The remaining acid is 75%.
    # The ratio [A-]/[HA] is therefore 0.25 / 0.75 = 1/3.
    # The volume term cancels out, so we can use the mole ratio directly.
    ratio_base_acid = 0.25 / 0.75
    calculated_ph_at_25_percent = pKa + math.log10(ratio_base_acid)

    # --- Calculation for pH at the Equivalence Point ---

    # 1. Calculate the initial moles of acetic acid before dilution.
    # This is the amount that will be titrated.
    V_acid_initial_L = V_acid_initial_cm3 / 1000.0
    moles_acid_initial = C_acid_initial_M * V_acid_initial_L

    # 2. Calculate the volume of NaOH needed to reach the equivalence point.
    # At equivalence, moles_base_added = moles_acid_initial.
    V_base_eq_L = moles_acid_initial / C_base_M

    # 3. Calculate the total volume of the solution at the equivalence point.
    # This is the volume of the diluted acid plus the volume of NaOH added.
    V_diluted_acid_L = (V_acid_initial_cm3 + V_water_cm3) / 1000.0
    V_total_at_eq_L = V_diluted_acid_L + V_base_eq_L

    # 4. Calculate the concentration of the conjugate base (acetate, A-) at the equivalence point.
    # All the initial acid has been converted to acetate.
    moles_acetate = moles_acid_initial
    C_acetate_M = moles_acetate / V_total_at_eq_L

    # 5. The acetate ion hydrolyzes water, producing OH-. We need Kb.
    # A- + H2O <=> HA + OH-
    Kb_acetate = Kw / Ka_acid

    # 6. Calculate the [OH-] concentration from the equilibrium.
    # Kb = [HA][OH-]/[A-] ≈ [OH-]^2 / C_acetate
    # We assume [OH-] is small compared to C_acetate.
    OH_concentration = math.sqrt(Kb_acetate * C_acetate_M)

    # 7. Calculate pOH and then the final pH.
    pOH = -math.log10(OH_concentration)
    calculated_ph_at_equivalence = 14.00 - pOH

    # --- Verification Step ---
    # Compare the calculated values with the expected values from Option C.
    # A tolerance is used for floating-point number comparison.
    tolerance = 0.015 # A slightly larger tolerance to account for rounding in the options

    is_ph_25_correct = math.isclose(calculated_ph_at_25_percent, expected_ph_at_25_percent, abs_tol=tolerance)
    is_ph_eq_correct = math.isclose(calculated_ph_at_equivalence, expected_ph_at_equivalence, abs_tol=tolerance)

    if is_ph_25_correct and is_ph_eq_correct:
        return "Correct"
    else:
        error_messages = []
        if not is_ph_25_correct:
            error_messages.append(
                f"Constraint for pH at 25% titration is not satisfied. "
                f"Calculated pH is {calculated_ph_at_25_percent:.2f}, but the answer provides {expected_ph_at_25_percent}."
            )
        if not is_ph_eq_correct:
            error_messages.append(
                f"Constraint for pH at the equivalence point is not satisfied. "
                f"Calculated pH is {calculated_ph_at_equivalence:.2f}, but the answer provides {expected_ph_at_equivalence}."
            )
        return "\n".join(error_messages)

# Run the check and print the result
result = check_answer()
print(result)