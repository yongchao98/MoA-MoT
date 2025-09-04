import math

def check_titration_ph_calculation():
    """
    This function verifies the pH calculations for the titration of diluted acetic acid with NaOH.
    It checks the pH at two specific points: 25% titration and the equivalence point.
    """
    # --- Problem Parameters ---
    V_acid_initial_cm3 = 20.00
    C_acid_initial = 0.05  # M (mol/L)
    V_water_cm3 = 20.00
    C_naoh = 0.1  # M (mol/L)
    Ka_acetic_acid = 1.85e-5
    Kw = 1.0e-14  # at 25 °C

    # --- Answer to Check (from option A) ---
    expected_ph_25_percent = 4.26
    expected_ph_equivalence = 8.52

    # --- Step 1: Calculate the concentration of acetic acid after dilution ---
    # Convert volumes to Liters
    V_acid_initial_L = V_acid_initial_cm3 / 1000.0
    V_water_L = V_water_cm3 / 1000.0

    # Calculate initial moles of acetic acid
    moles_acid = V_acid_initial_L * C_acid_initial

    # Calculate the final volume and concentration after dilution
    V_acid_diluted_L = V_acid_initial_L + V_water_L
    C_acid_diluted = moles_acid / V_acid_diluted_L

    # --- Step 2: Calculate pH at 25% titration ---
    # This is a buffer region, so the Henderson-Hasselbalch equation is used:
    # pH = pKa + log([A-]/[HA])
    # At 25% titration, 25% of HA is converted to A-, so 75% of HA remains.
    # The ratio [A-]/[HA] is 25/75 = 1/3.
    pKa = -math.log10(Ka_acetic_acid)
    ratio_25_percent = 0.25 / 0.75
    calculated_ph_25_percent = pKa + math.log10(ratio_25_percent)

    # --- Step 3: Calculate pH at the equivalence point ---
    # At the equivalence point, moles of NaOH added equals initial moles of acid.
    # All acetic acid is converted to its conjugate base, acetate (CH3COO-).
    
    # Volume of NaOH needed to reach equivalence point
    V_naoh_eq_L = moles_acid / C_naoh
    
    # Total volume of the solution at the equivalence point
    V_total_eq_L = V_acid_diluted_L + V_naoh_eq_L
    
    # Concentration of the acetate ion at the equivalence point
    C_acetate = moles_acid / V_total_eq_L
    
    # Acetate hydrolyzes in water: A- + H2O <=> HA + OH-
    # We need the base dissociation constant, Kb.
    Kb_acetate = Kw / Ka_acetic_acid
    
    # Use the equilibrium expression for the base: Kb = [HA][OH-]/[A-]
    # Let x = [OH-]. Then Kb = x^2 / (C_acetate - x).
    # We can approximate C_acetate - x ≈ C_acetate since Kb is small.
    # x^2 = Kb * C_acetate
    # x = [OH-]
    OH_concentration = math.sqrt(Kb_acetate * C_acetate)
    
    # Calculate pOH and then pH
    pOH = -math.log10(OH_concentration)
    calculated_ph_equivalence = 14.0 - pOH

    # --- Step 4: Compare calculated values with the expected answer ---
    # Use a small tolerance for floating-point comparisons.
    tolerance = 0.01

    is_25_percent_correct = math.isclose(calculated_ph_25_percent, expected_ph_25_percent, abs_tol=tolerance)
    is_equivalence_correct = math.isclose(calculated_ph_equivalence, expected_ph_equivalence, abs_tol=tolerance)

    if is_25_percent_correct and is_equivalence_correct:
        return "Correct"
    else:
        error_messages = []
        if not is_25_percent_correct:
            error_messages.append(
                f"The pH at 25% titration is incorrect. "
                f"Calculated value is {calculated_ph_25_percent:.2f}, but the answer states {expected_ph_25_percent}."
            )
        if not is_equivalence_correct:
            error_messages.append(
                f"The pH at the equivalence point is incorrect. "
                f"Calculated value is {calculated_ph_equivalence:.2f}, but the answer states {expected_ph_equivalence}."
            )
        return "\n".join(error_messages)

# Execute the check
result = check_titration_ph_calculation()
print(result)