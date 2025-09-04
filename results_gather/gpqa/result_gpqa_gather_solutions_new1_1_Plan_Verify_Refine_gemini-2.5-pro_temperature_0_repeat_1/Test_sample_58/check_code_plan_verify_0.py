import math

def check_ph_calculation():
    """
    This function checks the correctness of the pH calculation for the given mixture of acids and a base.
    It recalculates the pH based on the problem statement and compares it to the provided answer.
    """
    # --- Step 1: Define initial conditions from the question ---
    # CH3COOH (weak acid)
    vol_ch3cooh = 0.500  # L
    conc_ch3cooh = 0.1   # M

    # HCl (strong acid)
    vol_hcl = 0.400      # L
    conc_hcl = 0.2       # M

    # Ba(OH)2 (strong base)
    vol_baoh2 = 0.300    # L
    conc_baoh2 = 0.3     # M

    # The final answer provided is D, which corresponds to pH 12.62.
    expected_ph = 12.62

    # --- Step 2: Calculate initial moles of each species ---
    # Moles of weak acid
    moles_ch3cooh = vol_ch3cooh * conc_ch3cooh
    
    # Moles of H+ from strong acid (HCl)
    moles_h_from_hcl = vol_hcl * conc_hcl
    
    # Moles of OH- from strong base (Ba(OH)2)
    # A critical point: Ba(OH)2 is a dihydroxide base, providing 2 moles of OH- per mole.
    moles_oh_from_baoh2 = vol_baoh2 * conc_baoh2 * 2

    # --- Step 3: Perform neutralization calculation ---
    # The strong base (OH-) will react with all available acidic protons
    # from both the strong acid (HCl) and the weak acid (CH3COOH).
    total_moles_acid = moles_h_from_hcl + moles_ch3cooh
    total_moles_base = moles_oh_from_baoh2

    # Determine the excess species. The final pH will be determined by the excess strong acid or strong base.
    if total_moles_base > total_moles_acid:
        excess_moles_oh = total_moles_base - total_moles_acid
    else:
        # This case is not expected, but we handle it for completeness.
        excess_moles_h = total_moles_acid - total_moles_base
        # If there were excess acid, the calculation would be different, but the numbers show excess base.
        return f"Incorrect calculation: The result should be basic, but the calculation shows excess acid ({excess_moles_h:.4f} mol)."

    # --- Step 4: Calculate the final concentration of the excess species ---
    total_volume = vol_ch3cooh + vol_hcl + vol_baoh2
    final_oh_concentration = excess_moles_oh / total_volume

    # --- Step 5: Calculate the final pH ---
    if final_oh_concentration <= 0:
        return "Calculation error: Final OH- concentration is not a positive number."
        
    poh = -math.log10(final_oh_concentration)
    calculated_ph = 14 - poh

    # --- Step 6: Compare the calculated pH with the expected answer ---
    # Use a tolerance for floating-point comparison.
    tolerance = 0.01
    if abs(calculated_ph - expected_ph) < tolerance:
        return "Correct"
    else:
        reason = (
            f"The answer is incorrect. The calculated pH is {calculated_ph:.2f}, but the expected pH from the answer is {expected_ph}.\n"
            f"Here is the correct calculation breakdown:\n"
            f"1. Moles of CH3COOH = {moles_ch3cooh:.4f} mol.\n"
            f"2. Moles of H+ (from HCl) = {moles_h_from_hcl:.4f} mol.\n"
            f"3. Moles of OH- (from Ba(OH)2) = {moles_oh_from_baoh2:.4f} mol.\n"
            f"4. Total moles of acid (strong + weak) = {total_moles_acid:.4f} mol.\n"
            f"5. Total moles of base = {total_moles_base:.4f} mol.\n"
            f"6. Excess moles of OH- = {total_moles_base:.4f} - {total_moles_acid:.4f} = {excess_moles_oh:.4f} mol.\n"
            f"7. Total volume = {total_volume:.3f} L.\n"
            f"8. Final [OH-] = {excess_moles_oh:.4f} mol / {total_volume:.3f} L = {final_oh_concentration:.5f} M.\n"
            f"9. pOH = -log10({final_oh_concentration:.5f}) = {poh:.2f}.\n"
            f"10. pH = 14 - pOH = {calculated_ph:.2f}."
        )
        return reason

# Execute the check
result = check_ph_calculation()
print(result)