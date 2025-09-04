import math

def check_titration_answer():
    """
    This function calculates the pH at 25% titration and at the equivalence point
    for the given chemistry problem and compares it to the provided answer.
    """
    # --- Given values from the question ---
    V_HA_initial_cm3 = 20.00
    M_HA_initial = 0.05  # M (mol/L)
    V_water_cm3 = 20.00
    M_NaOH = 0.1  # M (mol/L)
    Ka = 1.85e-5
    Kw = 1.0e-14  # at 25 °C

    # --- Values from the provided answer (Option B) ---
    ans_pH_25 = 4.26
    ans_pH_eq = 8.52

    # --- Convert all volumes to Liters for consistency ---
    V_HA_initial_L = V_HA_initial_cm3 / 1000.0
    V_water_L = V_water_cm3 / 1000.0

    # --- Step 1: Calculate initial moles and concentration after dilution ---
    # Initial moles of acetic acid (HA)
    initial_moles_HA = M_HA_initial * V_HA_initial_L
    
    # Total volume after dilution with water
    V_diluted_L = V_HA_initial_L + V_water_L
    
    # Concentration of HA after dilution (this is the starting concentration for the titration)
    M_HA_diluted = initial_moles_HA / V_diluted_L

    # --- Step 2: Calculate pH at 25% titration ---
    # At 25% titration, 25% of the acetic acid has been neutralized to form its conjugate base, acetate (A-).
    # This creates a buffer solution where the ratio of [A-]/[HA] is 0.25 / (1 - 0.25) = 0.25 / 0.75 = 1/3.
    # We use the Henderson-Hasselbalch equation: pH = pKa + log([A-]/[HA])
    pKa = -math.log10(Ka)
    ratio_A_HA = 0.25 / 0.75
    calculated_pH_25 = pKa + math.log10(ratio_A_HA)

    # --- Step 3: Calculate pH at the equivalence point ---
    # At the equivalence point, all the acetic acid has been converted to acetate (A-).
    # Moles of NaOH added equals the initial moles of HA.
    moles_A_minus_at_eq = initial_moles_HA

    # Volume of NaOH needed to reach the equivalence point
    V_NaOH_eq_L = initial_moles_HA / M_NaOH
    
    # Total volume of the solution at the equivalence point
    V_total_eq_L = V_diluted_L + V_NaOH_eq_L
    
    # Concentration of the acetate ion (A-) at the equivalence point
    M_A_minus_at_eq = moles_A_minus_at_eq / V_total_eq_L

    # The acetate ion hydrolyzes in water: A- + H2O <=> HA + OH-
    # We need the base dissociation constant, Kb, for acetate.
    Kb = Kw / Ka
    
    # We can set up an ICE table for the hydrolysis, which gives Kb = [OH-]^2 / [A-].
    # We solve for [OH-], assuming [OH-] is much smaller than [A-].
    # [OH-] = sqrt(Kb * [A-])
    conc_OH = math.sqrt(Kb * M_A_minus_at_eq)
    
    # From [OH-], we calculate pOH and then pH.
    pOH = -math.log10(conc_OH)
    calculated_pH_eq = 14.0 - pOH

    # --- Step 4: Verify the answer ---
    # We check if the calculated values are close to the provided answer values, using a small tolerance.
    tolerance = 0.01
    is_pH_25_correct = abs(calculated_pH_25 - ans_pH_25) <= tolerance
    is_pH_eq_correct = abs(calculated_pH_eq - ans_pH_eq) <= tolerance

    if is_pH_25_correct and is_pH_eq_correct:
        return "Correct"
    else:
        reasons = []
        if not is_pH_25_correct:
            reasons.append(f"The pH at 25% titration is incorrect. Calculated value is {calculated_pH_25:.2f}, but the answer provided is {ans_pH_25}.")
        if not is_pH_eq_correct:
            reasons.append(f"The pH at the equivalence point is incorrect. Calculated value is {calculated_pH_eq:.2f}, but the answer provided is {ans_pH_eq}.")
        
        # Provide detailed calculation steps for clarity in case of an error.
        details = (
            f"\n\n--- Calculation Details ---\n"
            f"Initial moles of acetic acid: {initial_moles_HA:.5f} mol\n"
            f"Concentration after dilution: {M_HA_diluted:.3f} M\n"
            f"pKa of acetic acid: {pKa:.2f}\n"
            f"Calculated pH at 25% titration: {calculated_pH_25:.2f}\n"
            f"Volume of NaOH at equivalence: {V_NaOH_eq_L*1000:.2f} cm3\n"
            f"Total volume at equivalence: {V_total_eq_L*1000:.2f} cm3\n"
            f"Concentration of acetate [A-] at equivalence: {M_A_minus_at_eq:.4f} M\n"
            f"Kb of acetate: {Kb:.3e}\n"
            f"[OH-] at equivalence: {conc_OH:.3e} M\n"
            f"pOH at equivalence: {pOH:.2f}\n"
            f"Calculated pH at equivalence: {calculated_pH_eq:.2f}"
        )
        return "\n".join(reasons) + details

# To run the check, you would call the function:
# print(check_titration_answer())