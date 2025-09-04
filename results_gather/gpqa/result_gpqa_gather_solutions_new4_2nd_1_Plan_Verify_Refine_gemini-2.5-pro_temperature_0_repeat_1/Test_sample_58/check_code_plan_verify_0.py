import math

def check_answer():
    """
    This function checks the correctness of the pH calculation for the given mixture.
    """
    # 1. Define initial conditions from the question
    # CH3COOH (weak acid)
    vol_ch3cooh = 0.500  # L (500 mL)
    molarity_ch3cooh = 0.1  # M

    # HCl (strong acid)
    vol_hcl = 0.400  # L (400 mL)
    molarity_hcl = 0.2  # M

    # Ba(OH)2 (strong base)
    vol_baoh2 = 0.300  # L (300 mL)
    molarity_baoh2 = 0.3  # M

    # The final answer to check
    expected_ph = 12.62

    # 2. Calculate initial moles of acidic and basic species
    moles_ch3cooh = molarity_ch3cooh * vol_ch3cooh
    moles_h_plus_from_hcl = molarity_hcl * vol_hcl
    
    # Critical point: Ba(OH)2 provides 2 OH- ions per formula unit.
    moles_oh_minus_from_baoh2 = molarity_baoh2 * vol_baoh2 * 2

    # 3. Perform neutralization
    # The strong base will neutralize all available acidic protons from both strong and weak acids.
    total_moles_acid = moles_h_plus_from_hcl + moles_ch3cooh
    total_moles_base = moles_oh_minus_from_baoh2

    # 4. Determine the excess reactant and calculate its concentration
    if total_moles_base > total_moles_acid:
        # Base is in excess
        excess_moles_oh = total_moles_base - total_moles_acid
        
        # Calculate total volume
        total_volume = vol_ch3cooh + vol_hcl + vol_baoh2
        
        # Calculate final concentration of OH-
        final_concentration_oh = excess_moles_oh / total_volume
        
        # 5. Calculate pOH and then pH
        try:
            poh = -math.log10(final_concentration_oh)
            calculated_ph = 14 - poh
        except ValueError:
            return "Calculation error: Cannot take log of a non-positive number."

    elif total_moles_acid > total_moles_base:
        # Acid is in excess. This would require more complex buffer calculations if the
        # weak acid is left, but based on the numbers, this case shouldn't happen.
        # We will calculate the pH based on excess strong acid if that's the case.
        excess_moles_acid = total_moles_acid - total_moles_base
        total_volume = vol_ch3cooh + vol_hcl + vol_baoh2
        final_concentration_h = excess_moles_acid / total_volume
        try:
            calculated_ph = -math.log10(final_concentration_h)
        except ValueError:
            return "Calculation error: Cannot take log of a non-positive number."
    else:
        # Perfect neutralization
        calculated_ph = 7.0

    # 6. Compare the calculated pH with the expected answer
    tolerance = 0.01  # A small tolerance for floating point comparisons
    if abs(calculated_ph - expected_ph) < tolerance:
        return "Correct"
    else:
        # Provide a detailed reason for the discrepancy
        reason = (f"Incorrect. The calculated pH is {calculated_ph:.2f}, while the expected answer is {expected_ph}.\n"
                  f"Here's the breakdown:\n"
                  f"- Total moles of acid (H+ from HCl + CH3COOH): {total_moles_acid:.3f} mol\n"
                  f"- Total moles of base (OH- from Ba(OH)2): {total_moles_base:.3f} mol\n")
        if 'excess_moles_oh' in locals():
            reason += (f"- Excess moles of OH-: {excess_moles_oh:.3f} mol\n"
                       f"- Total Volume: {total_volume:.1f} L\n"
                       f"- Final [OH-]: {final_concentration_oh:.5f} M\n"
                       f"- Calculated pOH: {poh:.2f}\n"
                       f"- Calculated pH (14 - pOH): {calculated_ph:.2f}")
        return reason

# Run the check
print(check_answer())