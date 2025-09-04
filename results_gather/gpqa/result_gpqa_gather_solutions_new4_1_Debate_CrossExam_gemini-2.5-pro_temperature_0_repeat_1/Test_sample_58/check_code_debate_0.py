import math

def check_ph_calculation():
    """
    This function checks the correctness of the pH calculation for the given mixture.
    It recalculates the pH based on the problem's constraints and compares it
    to the provided answer.
    """
    
    # --- Problem Constraints & Provided Answer ---
    # Question: Determine the pH of the solution of 500 mL of 0.1 M CH3COOH, 
    # 400 mL of 0.2 M HCl, and 300 mL of 0.3 M Ba(OH)2.
    # Options: A) 12.62, B) 1.38, C) 8.68, D) 8.92
    # The provided answer to check is 'A'.
    
    options = {'A': 12.62, 'B': 1.38, 'C': 8.68, 'D': 8.92}
    provided_answer_key = 'A'
    expected_ph = options[provided_answer_key]

    # --- Step 1: Calculate initial moles of all species ---
    
    # Moles of CH3COOH (weak acid)
    vol_ch3cooh = 0.500  # L
    conc_ch3cooh = 0.1   # M
    moles_ch3cooh = vol_ch3cooh * conc_ch3cooh

    # Moles of H+ from HCl (strong acid)
    vol_hcl = 0.400  # L
    conc_hcl = 0.2   # M
    moles_h_strong = vol_hcl * conc_hcl

    # Moles of OH- from Ba(OH)2 (strong base)
    # Constraint check: Ba(OH)2 is a strong base that provides 2 OH- ions per formula unit.
    vol_baoh2 = 0.300  # L
    conc_baoh2 = 0.3   # M
    moles_oh_strong = vol_baoh2 * conc_baoh2 * 2

    # --- Step 2: Perform neutralization reaction ---
    
    # The strong base (OH-) will neutralize all available acid (strong first, then weak).
    # We can sum the total moles of acid to be neutralized.
    total_moles_acid = moles_h_strong + moles_ch3cooh
    
    # Compare total moles of base to total moles of acid
    if not math.isclose(total_moles_acid, 0.130):
        return f"Constraint check failed: Calculated total moles of acid is {total_moles_acid:.4f}, expected 0.130."
    if not math.isclose(moles_oh_strong, 0.180):
        return f"Constraint check failed: Calculated total moles of base is {moles_oh_strong:.4f}, expected 0.180."

    if moles_oh_strong > total_moles_acid:
        # Base is in excess, which determines the final pH.
        excess_moles_oh = moles_oh_strong - total_moles_acid
        
        # --- Step 3: Calculate final concentration of excess OH- ---
        
        # Constraint check: Total volume is the sum of the individual volumes.
        total_volume = vol_ch3cooh + vol_hcl + vol_baoh2
        if not math.isclose(total_volume, 1.2):
             return f"Constraint check failed: Calculated total volume is {total_volume:.4f} L, expected 1.2 L."

        final_conc_oh = excess_moles_oh / total_volume
        
        # --- Step 4: Calculate pOH and then pH ---
        
        pOH = -math.log10(final_conc_oh)
        calculated_ph = 14 - pOH
        
        # --- Step 5: Verify the answer ---
        
        # Check if the calculated pH matches the expected pH from the selected answer.
        if math.isclose(calculated_ph, expected_ph, rel_tol=1e-2):
            return "Correct"
        else:
            # If it doesn't match, explain why.
            reason = (
                f"The provided answer '{provided_answer_key}' ({expected_ph}) is incorrect.\n"
                f"The calculation shows that the final pH should be {calculated_ph:.2f}.\n"
                f"Here is the breakdown:\n"
                f"1. Total moles of acid (HCl + CH3COOH) = {total_moles_acid:.3f} mol.\n"
                f"2. Total moles of base (OH- from Ba(OH)2) = {moles_oh_strong:.3f} mol.\n"
                f"3. The strong base is in excess by {excess_moles_oh:.3f} mol.\n"
                f"4. The total volume is {total_volume:.1f} L.\n"
                f"5. The final concentration of excess [OH-] is {final_conc_oh:.5f} M.\n"
                f"6. This gives a pOH of {pOH:.2f}.\n"
                f"7. Therefore, the final pH is 14 - {pOH:.2f} = {calculated_ph:.2f}."
            )
            # Check if the calculated pOH matches another option, which is a common mistake.
            if math.isclose(pOH, options.get('B', -1), rel_tol=1e-2):
                reason += f"\nNote: The calculated pOH is {pOH:.2f}, which matches option B. The question asks for pH, not pOH."
            return reason
            
    elif total_moles_acid > moles_oh_strong:
        # This case is not expected for this problem but included for robustness.
        return "Incorrect: The calculation shows that acid is in excess, which would result in an acidic pH, contradicting the provided basic pH answer."
    else:
        # This case is for perfect neutralization.
        return "Incorrect: The calculation shows perfect neutralization (pH ~ 7.0), which contradicts the provided basic pH answer."

# Execute the check and print the result
result = check_ph_calculation()
print(result)