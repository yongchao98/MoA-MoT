import math

def check_ph_calculation():
    """
    This function checks the correctness of the pH calculation for the given mixture.
    It recalculates the pH from the initial conditions and compares it to the provided answer.
    """
    # --- Problem Constraints & Given Data ---
    # 500 mL of 0.1 M CH3COOH
    vol_ch3cooh = 0.500  # L
    conc_ch3cooh = 0.1   # M

    # 400 mL of 0.2 M HCl
    vol_hcl = 0.400      # L
    conc_hcl = 0.2       # M

    # 300 mL of 0.3 M Ba(OH)2
    vol_baoh2 = 0.300    # L
    conc_baoh2 = 0.3     # M

    # The provided answer is 'C', which corresponds to the value 12.62.
    expected_answer_value = 12.62

    # --- Step-by-step Calculation ---

    # 1. Calculate initial moles of all acidic and basic species.
    # Moles of weak acid (CH3COOH)
    moles_ch3cooh = vol_ch3cooh * conc_ch3cooh
    
    # Moles of H+ from strong acid (HCl)
    moles_h_from_hcl = vol_hcl * conc_hcl
    
    # Moles of OH- from strong base (Ba(OH)2).
    # It's crucial to remember that Ba(OH)2 is a dihydroxide base, providing 2 OH- ions.
    moles_oh_from_baoh2 = vol_baoh2 * conc_baoh2 * 2

    # 2. Determine total moles of acid and base to find the excess reactant.
    # The strong base will neutralize all available acid (both strong and weak).
    total_moles_acid = moles_h_from_hcl + moles_ch3cooh
    total_moles_base = moles_oh_from_baoh2

    # 3. Calculate the moles of the excess species.
    if total_moles_base <= total_moles_acid:
        return (f"Incorrect calculation: The base is not in excess as determined by the provided answer's logic. "
                f"Total moles of acid = {total_moles_acid:.3f}, Total moles of base = {total_moles_base:.3f}.")

    excess_moles_oh = total_moles_base - total_moles_acid

    # 4. Calculate the final concentration of the excess species.
    total_volume = vol_ch3cooh + vol_hcl + vol_baoh2
    final_conc_oh = excess_moles_oh / total_volume

    # 5. Calculate the final pH.
    # First, calculate pOH from the concentration of excess OH-.
    poh = -math.log10(final_conc_oh)
    # Then, calculate pH using the relation pH + pOH = 14.
    calculated_ph = 14 - poh

    # --- Verification ---
    # Compare the calculated pH with the expected answer value, allowing for a small tolerance.
    if math.isclose(calculated_ph, expected_answer_value, rel_tol=1e-2):
        return "Correct"
    else:
        reason = (
            f"The answer is incorrect.\n"
            f"The calculated pH is {calculated_ph:.2f}, but the provided answer 'C' corresponds to {expected_answer_value}.\n\n"
            f"Detailed calculation steps:\n"
            f"1. Total moles of acid = Moles(HCl) + Moles(CH3COOH) = {moles_h_from_hcl:.3f} + {moles_ch3cooh:.3f} = {total_moles_acid:.3f} mol.\n"
            f"2. Total moles of base (OH-) = Moles(Ba(OH)2) * 2 = ({vol_baoh2} L * {conc_baoh2} M) * 2 = {total_moles_base:.3f} mol.\n"
            f"3. Excess moles of OH- = {total_moles_base:.3f} - {total_moles_acid:.3f} = {excess_moles_oh:.3f} mol.\n"
            f"4. Total volume = {vol_ch3cooh} + {vol_hcl} + {vol_baoh2} = {total_volume:.1f} L.\n"
            f"5. Final [OH-] = {excess_moles_oh:.3f} mol / {total_volume:.1f} L = {final_conc_oh:.5f} M.\n"
            f"6. pOH = -log10({final_conc_oh:.5f}) = {poh:.2f}.\n"
            f"7. pH = 14 - pOH = 14 - {poh:.2f} = {calculated_ph:.2f}."
        )
        return reason

# Run the check and print the result
print(check_ph_calculation())