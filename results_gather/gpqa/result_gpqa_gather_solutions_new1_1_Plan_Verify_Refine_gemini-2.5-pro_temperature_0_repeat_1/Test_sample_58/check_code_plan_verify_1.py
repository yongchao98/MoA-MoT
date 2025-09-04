import math

def check_answer():
    """
    This function checks the correctness of the pH calculation for the given mixture.
    """
    # --- Define initial conditions from the question ---
    # CH3COOH (weak acid)
    vol_ch3cooh = 500 / 1000  # Convert mL to L
    conc_ch3cooh = 0.1        # M

    # HCl (strong acid)
    vol_hcl = 400 / 1000      # Convert mL to L
    conc_hcl = 0.2            # M

    # Ba(OH)2 (strong base)
    vol_baoh2 = 300 / 1000    # Convert mL to L
    conc_baoh2 = 0.3          # M

    # The final answer to check against
    expected_ph = 12.62

    # --- Step 1: Calculate initial moles of acidic and basic species ---
    moles_acid_weak = vol_ch3cooh * conc_ch3cooh
    moles_acid_strong = vol_hcl * conc_hcl
    
    # A critical constraint: Ba(OH)2 is a dihydroxide base, providing 2 OH- ions.
    moles_base_strong = vol_baoh2 * conc_baoh2 * 2

    # --- Step 2: Determine total moles for neutralization ---
    # A strong base will react with all available protons from both strong and weak acids.
    total_moles_acid = moles_acid_strong + moles_acid_weak
    total_moles_base = moles_base_strong

    # --- Step 3: Find the excess species after neutralization ---
    # The final pH will be determined by the excess strong acid or strong base.
    if total_moles_base > total_moles_acid:
        excess_moles_oh = total_moles_base - total_moles_acid
        
        # --- Step 4: Calculate final concentration and pH ---
        total_volume = vol_ch3cooh + vol_hcl + vol_baoh2
        final_conc_oh = excess_moles_oh / total_volume
        
        if final_conc_oh <= 0:
            return f"Calculation Error: Final hydroxide concentration is not positive ({final_conc_oh})."
            
        pOH = -math.log10(final_conc_oh)
        calculated_pH = 14 - pOH
    elif total_moles_acid > total_moles_base:
        # This case is not expected based on the numbers, but included for completeness.
        # The pH would be dominated by the remaining strong acid.
        moles_h_strong_remaining = moles_acid_strong - total_moles_base
        total_volume = vol_ch3cooh + vol_hcl + vol_baoh2
        final_conc_h = moles_h_strong_remaining / total_volume
        calculated_pH = -math.log10(final_conc_h)
    else:
        # Exact neutralization, resulting in a solution of the weak base CH3COO-.
        # This requires Ka/Kb for calculation and is not the case here.
        return "Calculation indicates exact neutralization, which is incorrect for the given values."

    # --- Step 5: Verify the result ---
    tolerance = 0.01
    if abs(calculated_pH - expected_ph) < tolerance:
        return "Correct"
    else:
        return (
            f"Incorrect. The calculated pH is {calculated_pH:.2f}, which does not match the expected pH of {expected_ph}.\n"
            f"Calculation breakdown:\n"
            f" - Total moles of acid (HCl + CH3COOH): {total_moles_acid:.4f} mol\n"
            f" - Total moles of base (OH- from Ba(OH)2): {total_mles_base:.4f} mol\n"
            f" - Excess moles of OH-: {excess_moles_oh:.4f} mol\n"
            f" - Total Volume: {total_volume:.3f} L\n"
            f" - Final [OH-]: {final_conc_oh:.5f} M\n"
            f" - pOH: {-math.log10(final_conc_oh):.2f}\n"
            f" - Final pH: {calculated_pH:.2f}"
        )

# Run the check
result = check_answer()
print(result)