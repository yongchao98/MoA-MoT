import math

def check_ph_calculation():
    """
    This function checks the correctness of the pH calculation for the given mixture.
    The problem is to find the pH of a solution made by mixing:
    - 500 mL of 0.1 M CH3COOH (weak acid)
    - 400 mL of 0.2 M HCl (strong acid)
    - 300 mL of 0.3 M Ba(OH)2 (strong base)
    The expected correct answer from the options is pH = 12.62 (Option C).
    This code will perform the calculation and verify if that answer is correct.
    """
    
    # --- Step 1: Calculate initial moles of all acidic and basic species ---
    
    # Moles of CH3COOH (weak acid)
    vol_ch3cooh = 0.500  # L
    conc_ch3cooh = 0.1   # M
    moles_ch3cooh = vol_ch3cooh * conc_ch3cooh

    # Moles of H+ from HCl (strong acid)
    vol_hcl = 0.400  # L
    conc_hcl = 0.2   # M
    moles_h_strong = vol_hcl * conc_hcl

    # Moles of OH- from Ba(OH)2 (strong base)
    vol_baoh2 = 0.300  # L
    conc_baoh2 = 0.3   # M
    # Ba(OH)2 provides 2 OH- ions per formula unit
    moles_oh_strong = 2 * (vol_baoh2 * conc_baoh2)

    # --- Step 2: Simulate neutralization reactions ---
    # The strong base (OH-) reacts completely with the strong acid (H+) first.
    moles_oh_after_strong_acid_neut = moles_oh_strong - moles_h_strong

    # Now, the remaining strong base (OH-) reacts with the weak acid (CH3COOH).
    moles_oh_final = moles_oh_after_strong_acid_neut - moles_ch3cooh

    # --- Step 3: Calculate the final pH ---
    # The final solution contains an excess of a strong base (OH-).
    # The pH will be determined by the concentration of this excess OH-.
    
    # Calculate total volume
    total_volume = vol_ch3cooh + vol_hcl + vol_baoh2

    # Check if there is excess base. If moles_oh_final is positive, there is.
    if moles_oh_final > 0:
        # Calculate the concentration of the excess OH-
        final_conc_oh = moles_oh_final / total_volume
        
        # Calculate pOH and then pH
        pOH = -math.log10(final_conc_oh)
        calculated_pH = 14.0 - pOH
    else:
        # This block would handle cases with excess acid or a buffer,
        # but for this specific problem, we have excess strong base.
        # We can return an error if this unexpected state is reached.
        return "Calculation error: The final solution was not basic as expected."

    # --- Step 4: Check against the provided answer ---
    # The correct option is C, which is 12.62.
    expected_pH = 12.62
    
    # We use a small tolerance for floating point comparison
    tolerance = 0.01
    
    if abs(calculated_pH - expected_pH) < tolerance:
        return "Correct"
    else:
        return (f"Incorrect. The calculated pH is {calculated_pH:.2f}, which does not match the expected answer of {expected_pH}. "
                f"The calculation shows that after all neutralizations, there are {moles_oh_final:.2f} moles of excess OH- in a total volume of {total_volume:.1f} L, "
                f"leading to a final [OH-] of {final_conc_oh:.4f} M and a pOH of {pOH:.2f}.")

# To check the answer, we run the function.
# The function will return "Correct" if the calculation confirms that 12.62 is the right pH.
result = check_ph_calculation()
# print(result) # This would print "Correct"