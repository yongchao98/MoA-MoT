import math

def check_ph_calculation():
    """
    This function calculates the pH of the given solution mixture and checks if the provided answer is correct.
    
    The mixture consists of:
    - 500 mL of 0.1 M CH3COOH (weak acid)
    - 400 mL of 0.2 M HCl (strong acid)
    - 300 mL of 0.3 M Ba(OH)2 (strong base)
    
    The provided answer is B, which corresponds to a pH of 12.62.
    """
    
    # --- Step 1: Define initial conditions and calculate moles ---
    
    # Acetic Acid (weak acid)
    vol_ch3cooh = 0.500  # L
    conc_ch3cooh = 0.1   # M
    moles_ch3cooh = vol_ch3cooh * conc_ch3cooh
    
    # Hydrochloric Acid (strong acid)
    vol_hcl = 0.400  # L
    conc_hcl = 0.2   # M
    # HCl is a strong acid, so moles H+ = moles HCl
    moles_h_strong = vol_hcl * conc_hcl
    
    # Barium Hydroxide (strong base)
    vol_baoh2 = 0.300  # L
    conc_baoh2 = 0.3   # M
    moles_baoh2 = vol_baoh2 * conc_baoh2
    # Ba(OH)2 provides 2 OH- ions per formula unit
    moles_oh_strong = moles_baoh2 * 2
    
    # --- Step 2: Perform neutralization reaction ---
    
    # The strong base (OH-) will neutralize all available acidic protons,
    # both from the strong acid (HCl) and the weak acid (CH3COOH).
    total_moles_acid = moles_h_strong + moles_ch3cooh
    
    # Compare total moles of acid and base to find the excess reactant
    if moles_oh_strong > total_moles_acid:
        # Strong base is in excess
        excess_moles_oh = moles_oh_strong - total_moles_acid
        
        # --- Step 3: Calculate final concentration and pH ---
        
        # Calculate total volume of the final solution
        total_volume = vol_ch3cooh + vol_hcl + vol_baoh2
        
        # Calculate the concentration of the excess OH- ions
        final_conc_oh = excess_moles_oh / total_volume
        
        # Calculate pOH
        pOH = -math.log10(final_conc_oh)
        
        # Calculate pH
        pH = 14 - pOH
        
    elif total_moles_acid > moles_oh_strong:
        # Acid is in excess. This case is more complex as it could result
        # in a buffer solution plus excess strong acid.
        # However, based on the problem's options, excess base is expected.
        # We can calculate the remaining acid species for completeness.
        excess_moles_acid = total_moles_acid - moles_oh_strong
        # This case is not expected, so we'll return an error if it occurs.
        return f"Incorrect: Calculation shows excess acid, which does not match the provided answer's logic. Excess acid moles: {excess_moles_acid:.4f}"
        
    else:
        # Perfect neutralization. The pH would be determined by the hydrolysis
        # of the acetate ion (CH3COO-), resulting in a slightly basic solution.
        # This is also not expected.
        return "Incorrect: Calculation shows perfect neutralization, which does not match the provided answer's logic."

    # --- Step 4: Check the correctness of the answer ---
    
    # The provided answer is B, which is 12.62
    expected_ph = 12.62
    
    # Check if the calculated pH is close to the expected pH (within a tolerance)
    if abs(pH - expected_ph) < 0.01:
        return "Correct"
    else:
        return (f"Incorrect: The calculated pH is {pH:.2f}, but the expected answer is {expected_ph}. "
                f"The calculation steps are as follows:\n"
                f"1. Total moles of acid (H+ from HCl + CH3COOH) = {moles_h_strong:.3f} + {moles_ch3cooh:.3f} = {total_moles_acid:.3f} mol.\n"
                f"2. Total moles of base (OH- from Ba(OH)2) = {moles_baoh2:.3f} * 2 = {moles_oh_strong:.3f} mol.\n"
                f"3. Excess moles of OH- = {moles_oh_strong:.3f} - {total_moles_acid:.3f} = {excess_moles_oh:.3f} mol.\n"
                f"4. Total volume = {total_volume:.1f} L.\n"
                f"5. Final [OH-] = {excess_moles_oh:.3f} / {total_volume:.1f} = {final_conc_oh:.5f} M.\n"
                f"6. pOH = -log10({final_conc_oh:.5f}) = {pOH:.2f}.\n"
                f"7. pH = 14 - pOH = {pH:.2f}.")

# Run the check
result = check_ph_calculation()
print(result)