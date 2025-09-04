import math

def check_ph_correctness():
    """
    This function checks the correctness of the given answer for the pH calculation problem.
    """
    # --- Input values from the question ---
    # Acetic Acid (weak acid)
    vol_ch3cooh = 0.500  # L
    conc_ch3cooh = 0.1   # M

    # Hydrochloric Acid (strong acid)
    vol_hcl = 0.400      # L
    conc_hcl = 0.2       # M

    # Barium Hydroxide (strong base)
    vol_baoh2 = 0.300    # L
    conc_baoh2 = 0.3     # M

    # The answer provided by the LLM to be checked
    llm_answer_ph = 12.62
    llm_answer_option = 'B'

    # --- Step 1: Calculate initial moles of reactive species ---
    moles_h_strong = vol_hcl * conc_hcl
    moles_weak_acid = vol_ch3cooh * conc_ch3cooh
    # Ba(OH)2 is a strong base that provides 2 OH- ions per formula unit
    moles_oh_strong = vol_baoh2 * conc_baoh2 * 2

    # --- Step 2: Perform neutralization ---
    # The strong base will neutralize all available acid (strong first, then weak)
    total_moles_acid = moles_h_strong + moles_weak_acid
    
    if moles_oh_strong > total_moles_acid:
        # Excess strong base determines the pH
        excess_moles_oh = moles_oh_strong - total_moles_acid
        total_volume = vol_ch3cooh + vol_hcl + vol_baoh2
        final_oh_conc = excess_moles_oh / total_volume
        poh = -math.log10(final_oh_conc)
        calculated_ph = 14 - poh
    elif total_moles_acid > moles_oh_strong:
        # Excess acid determines the pH. We must check if it's strong or weak acid remaining.
        moles_oh_after_strong_acid_reaction = moles_oh_strong - moles_h_strong
        if moles_oh_after_strong_acid_reaction < 0: # Excess strong acid
            excess_moles_h = abs(moles_oh_after_strong_acid_reaction)
            total_volume = vol_ch3cooh + vol_hcl + vol_baoh2
            final_h_conc = excess_moles_h / total_volume
            calculated_ph = -math.log10(final_h_conc)
        else: # Buffer solution or weak acid solution - more complex calculation needed
            calculated_ph = -1 # Placeholder for more complex cases
    else:
        # Exact neutralization leading to a salt solution
        calculated_ph = 7 # Placeholder, would depend on salt hydrolysis

    # --- Step 3: Compare calculated result with the LLM's answer ---
    # Use a small tolerance for floating point comparison
    if math.isclose(calculated_ph, llm_answer_ph, rel_tol=1e-2):
        return "Correct"
    else:
        return (f"Incorrect. The provided answer is {llm_answer_ph} (Option {llm_answer_option}), but the calculated pH is {calculated_ph:.2f}.\n"
                f"Reasoning:\n"
                f"1. Moles of H+ (from HCl) = {moles_h_strong:.3f} mol.\n"
                f"2. Moles of CH3COOH = {moles_weak_acid:.3f} mol.\n"
                f"3. Total moles of acid = {total_moles_acid:.3f} mol.\n"
                f"4. Moles of OH- (from Ba(OH)2) = {moles_oh_strong:.3f} mol.\n"
                f"5. There is an excess of strong base (OH-). Excess moles = {moles_oh_strong:.3f} - {total_moles_acid:.3f} = {excess_moles_oh:.3f} mol.\n"
                f"6. Total volume = {total_volume:.1f} L.\n"
                f"7. Final [OH-] = {excess_moles_oh:.3f} / {total_volume:.1f} = {final_oh_conc:.4f} M.\n"
                f"8. pOH = -log10({final_oh_conc:.4f}) = {-math.log10(final_oh_conc):.2f}.\n"
                f"9. pH = 14 - pOH = {14 - (-math.log10(final_oh_conc)):.2f}.")

# Execute the check
result = check_ph_correctness()
print(result)