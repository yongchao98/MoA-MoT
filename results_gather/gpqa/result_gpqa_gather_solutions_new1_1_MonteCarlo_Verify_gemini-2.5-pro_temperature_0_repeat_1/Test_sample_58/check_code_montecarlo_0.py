import math

def check_answer():
    """
    This function checks the correctness of the provided answer for the chemistry pH problem.
    """
    # --- Problem Data ---
    # CH3COOH (weak acid)
    vol_ch3cooh_L = 500 / 1000
    molarity_ch3cooh = 0.1

    # HCl (strong acid)
    vol_hcl_L = 400 / 1000
    molarity_hcl = 0.2

    # Ba(OH)2 (strong base)
    vol_baoh2_L = 300 / 1000
    molarity_baoh2 = 0.3

    # --- Provided Answer ---
    # The final answer from the analysis is B, which corresponds to pH = 12.62.
    expected_ph = 12.62
    expected_option = 'B'

    # --- Calculation ---
    # Step 1: Calculate initial moles of all acidic and basic species.
    # For neutralization with a strong base, we consider all acidic protons from both weak and strong acids.
    moles_acid_ch3cooh = vol_ch3cooh_L * molarity_ch3cooh
    moles_acid_hcl = vol_hcl_L * molarity_hcl
    total_moles_acid = moles_acid_ch3cooh + moles_acid_hcl

    # Ba(OH)2 is a strong base that provides 2 moles of OH- for every 1 mole of Ba(OH)2.
    moles_base_baoh2 = vol_baoh2_L * molarity_baoh2
    total_moles_base = moles_base_baoh2 * 2

    # Step 2: Determine the excess species after neutralization.
    if total_moles_base > total_moles_acid:
        excess_moles_oh = total_moles_base - total_moles_acid
        
        # Step 3: Calculate the total volume of the final solution.
        total_volume_L = vol_ch3cooh_L + vol_hcl_L + vol_baoh2_L

        # Step 4: Calculate the concentration of the excess OH- ions.
        final_conc_oh = excess_moles_oh / total_volume_L

        # Step 5: Calculate pOH and then pH.
        pOH = -math.log10(final_conc_oh)
        calculated_ph = 14 - pOH
    elif total_moles_acid > total_moles_base:
        # This case is not expected based on the options, but we handle it for completeness.
        # The strong base neutralizes the strong acid first.
        moles_hcl_remaining = moles_acid_hcl - total_moles_base
        if moles_hcl_remaining > 0: # Strong acid is still in excess
            total_volume_L = vol_ch3cooh_L + vol_hcl_L + vol_baoh2_L
            final_conc_h = moles_hcl_remaining / total_volume_L
            calculated_ph = -math.log10(final_conc_h)
        else: # All strong acid is gone, some weak acid is neutralized -> buffer
            # This scenario would lead to a buffer calculation, which is more complex
            # and doesn't match the simple numeric options.
            return "Calculation error: The reaction results in a buffer solution, which does not align with the provided options."
    else:
        # Perfect neutralization. The solution contains the conjugate base of a weak acid.
        # pH would be > 7 but requires a Kb calculation. This is not the expected path.
        return "Calculation error: Perfect neutralization occurred, which does not align with the provided options."

    # --- Verification ---
    # Check if the calculated pH matches the expected pH within a small tolerance.
    tolerance = 0.01
    if abs(calculated_ph - expected_ph) < tolerance:
        return "Correct"
    else:
        error_message = (
            f"The provided answer is {expected_option} (pH={expected_ph}), but the calculation yields a different result.\n"
            f"Step-by-step calculation:\n"
            f"1. Total moles of acid (H+ from HCl + CH3COOH) = {moles_acid_hcl:.4f} + {moles_acid_ch3cooh:.4f} = {total_moles_acid:.4f} mol\n"
            f"2. Total moles of base (OH- from Ba(OH)2) = {moles_base_baoh2:.4f} mol Ba(OH)2 * 2 = {total_moles_base:.4f} mol\n"
            f"3. Excess species: Base is in excess by {excess_moles_oh:.4f} mol.\n"
            f"4. Total volume = {total_volume_L:.3f} L\n"
            f"5. Final [OH-] concentration = {excess_moles_oh:.4f} mol / {total_volume_L:.3f} L = {final_conc_oh:.5f} M\n"
            f"6. Calculated pOH = -log10({final_conc_oh:.5f}) = {pOH:.2f}\n"
            f"7. Calculated pH = 14 - {pOH:.2f} = {calculated_ph:.2f}\n"
            f"The calculated pH of {calculated_ph:.2f} does not match the expected pH of {expected_ph}."
        )
        return error_message

# Run the check and print the result
print(check_answer())