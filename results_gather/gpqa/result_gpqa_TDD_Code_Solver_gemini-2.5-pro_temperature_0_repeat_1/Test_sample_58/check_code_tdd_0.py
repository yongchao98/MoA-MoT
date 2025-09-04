import math

def check_ph_calculation():
    """
    This function checks the correctness of the pH calculation for the given mixture.
    The mixture consists of:
    - 500 mL of 0.1 M CH3COOH (a weak acid)
    - 400 mL of 0.2 M HCl (a strong acid)
    - 300 mL of 0.3 M Ba(OH)2 (a strong base)
    The provided answer is C, which corresponds to a pH of 12.62.
    """
    # --- Problem Data and Provided Answer ---
    # Acetic Acid (Weak Acid)
    vol_ch3cooh = 0.500  # L
    conc_ch3cooh = 0.1   # M
    # Hydrochloric Acid (Strong Acid)
    vol_hcl = 0.400      # L
    conc_hcl = 0.2       # M
    # Barium Hydroxide (Strong Base)
    vol_baoh2 = 0.300    # L
    conc_baoh2 = 0.3     # M

    # The answer to check, corresponding to option C
    expected_ph = 12.62

    # --- Step 1: Calculate initial moles of acid and base ---
    # The strong base will react with all available protons, from both strong and weak acids.
    moles_acid_total = (vol_ch3cooh * conc_ch3cooh) + (vol_hcl * conc_hcl)
    
    # Ba(OH)2 is a strong base that dissociates to yield two OH- ions per formula unit.
    moles_base_total = (vol_baoh2 * conc_baoh2) * 2

    # --- Step 2: Determine excess moles after neutralization ---
    # The reaction between H+ and OH- goes to completion.
    if moles_base_total > moles_acid_total:
        excess_moles_oh = moles_base_total - moles_acid_total
    else:
        # This case would result in an acidic solution, which is not what happens here.
        # For completeness, we calculate the excess acid.
        excess_moles_h = moles_acid_total - moles_base_total
        # If this path were taken, the answer would be incorrect.
        # We can construct an error message here, but the main check is at the end.

    # --- Step 3: Calculate final pH ---
    # The final pH is determined by the species in excess.
    total_volume = vol_ch3cooh + vol_hcl + vol_baoh2
    
    if 'excess_moles_oh' in locals():
        # The solution is basic due to excess strong base.
        final_oh_concentration = excess_moles_oh / total_volume
        poh = -math.log10(final_oh_concentration)
        calculated_ph = 14.0 - poh
    elif 'excess_moles_h' in locals():
        # The solution would be acidic.
        final_h_concentration = excess_moles_h / total_volume
        calculated_ph = -math.log10(final_h_concentration)
    else: # moles_base_total == moles_acid_total
        # This case would require calculating pH from the salt of a weak acid.
        # Not applicable for this problem.
        return "Calculation error: Perfect neutralization is not expected."

    # --- Step 4: Verify the answer ---
    # Compare the calculated pH with the expected pH, allowing for a small rounding tolerance.
    if abs(calculated_ph - expected_ph) < 0.01:
        return "Correct"
    else:
        # If the answer is incorrect, provide a detailed breakdown of the correct calculation.
        reason = (
            f"The provided answer pH of {expected_ph} is incorrect. The calculated pH is {calculated_ph:.2f}.\n\n"
            f"Calculation Breakdown:\n"
            f"1. Total Moles of Acid (H+):\n"
            f"   - Moles from CH3COOH = 0.500 L * 0.1 M = {vol_ch3cooh * conc_ch3cooh:.3f} mol\n"
            f"   - Moles from HCl = 0.400 L * 0.2 M = {vol_hcl * conc_hcl:.3f} mol\n"
            f"   - Total Acid = {vol_ch3cooh * conc_ch3cooh:.3f} + {vol_hcl * conc_hcl:.3f} = {moles_acid_total:.3f} mol\n\n"
            f"2. Total Moles of Base (OH-):\n"
            f"   - Moles from Ba(OH)2 = 0.300 L * 0.3 M * 2 = {moles_base_total:.3f} mol\n\n"
            f"3. Net Moles after Neutralization:\n"
            f"   - Since Moles OH- ({moles_base_total:.3f}) > Moles H+ ({moles_acid_total:.3f}), there is excess OH-.\n"
            f"   - Excess Moles OH- = {moles_base_total:.3f} - {moles_acid_total:.3f} = {excess_moles_oh:.3f} mol\n\n"
            f"4. Final pH Calculation:\n"
            f"   - Total Volume = 0.500 L + 0.400 L + 0.300 L = {total_volume:.1f} L\n"
            f"   - Final [OH-] = {excess_moles_oh:.3f} mol / {total_volume:.1f} L = {final_oh_concentration:.5f} M\n"
            f"   - pOH = -log10({final_oh_concentration:.5f}) = {-math.log10(final_oh_concentration):.2f}\n"
            f"   - pH = 14 - pOH = {calculated_ph:.2f}"
        )
        return reason

# Execute the check and print the result
result = check_ph_calculation()
print(result)