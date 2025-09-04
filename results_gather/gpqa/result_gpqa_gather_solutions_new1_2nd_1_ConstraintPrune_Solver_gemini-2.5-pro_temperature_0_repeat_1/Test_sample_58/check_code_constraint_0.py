import math

def check_ph_calculation():
    """
    This function checks the correctness of the pH calculation for the given mixture of solutions.
    
    Question: Determine the pH of the solution of 500 mL of 0.1 M CH3COOH, 
              400 mL of 0.2 M HCl, and 300 mL of 0.3 M Ba(OH)2.
    Options: A) 8.92, B) 1.38, C) 8.68, D) 12.62
    The provided answer is D, which corresponds to a pH of 12.62.
    """
    
    # --- Step 1: Define initial conditions and calculate moles ---
    
    # CH3COOH (weak acid)
    vol_ch3cooh = 0.500  # L
    conc_ch3cooh = 0.1   # M
    moles_ch3cooh = vol_ch3cooh * conc_ch3cooh
    
    # HCl (strong acid)
    vol_hcl = 0.400  # L
    conc_hcl = 0.2   # M
    moles_hcl = vol_hcl * conc_hcl
    
    # Ba(OH)2 (strong base)
    vol_baoh2 = 0.300  # L
    conc_baoh2 = 0.3   # M
    moles_baoh2 = vol_baoh2 * conc_baoh2
    
    # --- Step 2: Calculate total moles of acid and base ---
    
    # Total moles of acid (protons available for neutralization)
    # The strong base will neutralize both the strong acid and the weak acid.
    total_moles_acid = moles_hcl + moles_ch3cooh
    
    # Total moles of base (hydroxide ions)
    # Ba(OH)2 is a dihydroxide base, so it provides 2 moles of OH- per mole of Ba(OH)2.
    moles_oh = moles_baoh2 * 2
    
    # --- Step 3: Perform neutralization and find the excess reactant ---
    
    if moles_oh > total_moles_acid:
        # Base is in excess
        excess_moles_oh = moles_oh - total_moles_acid
        
        # --- Step 4: Calculate final concentration of the excess ion ---
        total_volume = vol_ch3cooh + vol_hcl + vol_baoh2
        final_conc_oh = excess_moles_oh / total_volume
        
        # --- Step 5: Calculate pOH and then pH ---
        try:
            poh = -math.log10(final_conc_oh)
            calculated_ph = 14 - poh
        except ValueError:
            return "Error: Cannot calculate log of a non-positive number."
            
    elif total_moles_acid > moles_oh:
        # Acid is in excess
        excess_moles_h = total_moles_acid - moles_oh
        
        # --- Step 4: Calculate final concentration of the excess ion ---
        total_volume = vol_ch3cooh + vol_hcl + vol_baoh2
        final_conc_h = excess_moles_h / total_volume
        
        # --- Step 5: Calculate pH ---
        try:
            calculated_ph = -math.log10(final_conc_h)
        except ValueError:
            return "Error: Cannot calculate log of a non-positive number."
            
    else:
        # Neutralization is complete. This case is more complex, but based on the numbers, it's unlikely.
        # The solution would contain the conjugate base CH3COO-, making it slightly basic.
        # However, the problem setup points to an excess of a strong species.
        calculated_ph = 7.0 # Placeholder for perfect neutralization of strong acid/base
        # A more detailed calculation would be needed for the weak base solution.

    # --- Step 6: Compare the calculated pH with the provided answer ---
    
    expected_ph = 12.62 # This corresponds to option D
    
    # Check if the calculated pH is very close to the expected pH
    if abs(calculated_ph - expected_ph) < 0.01:
        return "Correct"
    else:
        # Provide a detailed reason for the discrepancy
        reason = f"The answer is incorrect.\n"
        reason += f"Calculation details:\n"
        reason += f"  - Total moles of acid (H+): {total_moles_acid:.4f} mol\n"
        reason += f"  - Total moles of base (OH-): {moles_oh:.4f} mol\n"
        if 'excess_moles_oh' in locals():
            reason += f"  - Excess moles of OH-: {excess_moles_oh:.4f} mol\n"
            reason += f"  - Total volume: {total_volume:.1f} L\n"
            reason += f"  - Final [OH-]: {final_conc_oh:.5f} M\n"
            reason += f"  - Calculated pOH: {-math.log10(final_conc_oh):.2f}\n"
        elif 'excess_moles_h' in locals():
            reason += f"  - Excess moles of H+: {excess_moles_h:.4f} mol\n"
            reason += f"  - Total volume: {total_volume:.1f} L\n"
            reason += f"  - Final [H+]: {final_conc_h:.5f} M\n"
        
        reason += f"The calculated pH is {calculated_ph:.2f}, which does not match the expected pH of {expected_ph}.\n"
        
        # Check for common errors
        if 'poh' in locals() and abs(poh - 1.38) < 0.01:
             reason += "The calculated pOH is 1.38. The error in some answers might be confusing pOH (option B) with pH."
        
        return reason

# Run the check and print the result
print(check_ph_calculation())