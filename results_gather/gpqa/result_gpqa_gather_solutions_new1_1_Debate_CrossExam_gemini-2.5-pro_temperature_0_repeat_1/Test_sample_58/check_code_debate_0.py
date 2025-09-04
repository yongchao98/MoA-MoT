import math

def check_ph_calculation():
    """
    This function checks the correctness of the pH calculation for the given chemical mixture.
    
    The problem involves mixing:
    - 500 mL of 0.1 M CH3COOH (weak acid)
    - 400 mL of 0.2 M HCl (strong acid)
    - 300 mL of 0.3 M Ba(OH)2 (strong base)
    
    The expected answer is D, which corresponds to a pH of 12.62.
    """
    
    # --- Step 1: Define initial conditions and calculate moles ---
    
    # Acetic Acid (CH3COOH) - weak acid
    vol_ch3cooh = 0.500  # L
    conc_ch3cooh = 0.1   # M
    moles_ch3cooh = vol_ch3cooh * conc_ch3cooh
    
    # Hydrochloric Acid (HCl) - strong acid
    vol_hcl = 0.400  # L
    conc_hcl = 0.2   # M
    moles_h_from_hcl = vol_hcl * conc_hcl
    
    # Barium Hydroxide (Ba(OH)2) - strong base
    vol_baoh2 = 0.300  # L
    conc_baoh2 = 0.3   # M
    moles_baoh2 = vol_baoh2 * conc_baoh2
    # Crucial point: Ba(OH)2 is a dihydroxide base, so it provides 2 moles of OH- per mole of Ba(OH)2.
    moles_oh_from_baoh2 = moles_baoh2 * 2
    
    # --- Step 2: Perform neutralization reaction ---
    
    # The strong base (OH-) will neutralize all available acidic protons, both from the strong acid (HCl)
    # and the weak acid (CH3COOH).
    total_moles_acid = moles_h_from_hcl + moles_ch3cooh
    total_moles_base = moles_oh_from_baoh2
    
    # Check if the final solution is acidic or basic
    if total_moles_base > total_moles_acid:
        # Base is in excess
        excess_moles_oh = total_moles_base - total_moles_acid
        
        # --- Step 3: Calculate final concentration and pH ---
        
        # Calculate total volume
        total_volume = vol_ch3cooh + vol_hcl + vol_baoh2
        
        # Calculate final concentration of excess OH-
        final_conc_oh = excess_moles_oh / total_volume
        
        # Calculate pOH
        try:
            poh = -math.log10(final_conc_oh)
        except ValueError:
            return "Error: Cannot calculate log of a non-positive concentration."
            
        # Calculate pH
        ph = 14 - poh
        
        # --- Step 4: Verify the answer ---
        
        # The expected pH from option D is 12.62
        expected_ph = 12.62
        
        if math.isclose(ph, expected_ph, rel_tol=1e-2):
            return "Correct"
        else:
            return (f"Incorrect. The calculation shows an excess of {excess_moles_oh:.4f} moles of OH- in a total volume of {total_volume:.1f} L, "
                    f"leading to [OH-] = {final_conc_oh:.5f} M. This gives a pOH of {poh:.2f} and a final pH of {ph:.2f}. "
                    f"This matches the value of 12.62 for option D, but the provided answer might have a different reasoning or final choice.")

    elif total_moles_acid > total_moles_base:
        # Acid is in excess. This would be a more complex calculation involving either a buffer or a strong acid solution.
        # However, based on the problem's options, this outcome is not expected.
        return (f"Incorrect. The calculation shows that the acid is in excess, which would result in an acidic pH. "
                f"Total moles of acid = {total_moles_acid:.4f}, Total moles of base = {total_moles_base:.4f}. "
                f"This contradicts the provided answer which is strongly basic.")
    else:
        # Neutral solution (pH=7), also not expected.
        return "Incorrect. The calculation shows that the moles of acid and base are exactly equal, which would result in a neutral solution (or one determined by the weak conjugate base), not a pH of 12.62."

# Run the check
result = check_ph_calculation()
print(result)