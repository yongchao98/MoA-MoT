import math

def check_correctness():
    """
    This function checks the correctness of the given LLM's answer for a chemistry pH calculation problem.
    """
    # --- Problem Data ---
    # Solution 1: Weak Acid
    vol_ch3cooh = 0.500  # L (500 mL)
    conc_ch3cooh = 0.1   # M
    # Solution 2: Strong Acid
    vol_hcl = 0.400      # L (400 mL)
    conc_hcl = 0.2       # M
    # Solution 3: Strong Base
    vol_baoh2 = 0.300    # L (300 mL)
    conc_baoh2 = 0.3     # M

    # --- LLM's Answer ---
    llm_answer_ph = 12.62

    # --- Calculation ---

    # Step 1: Calculate initial moles of all reactants.
    moles_ch3cooh = vol_ch3cooh * conc_ch3cooh
    moles_hcl = vol_hcl * conc_hcl
    moles_baoh2 = vol_baoh2 * conc_baoh2

    # Step 2: Calculate initial moles of H+ from the strong acid and OH- from the strong base.
    # HCl is a strong acid, so moles H+ = moles HCl.
    moles_h_strong = moles_hcl
    # Ba(OH)2 is a strong base that provides 2 moles of OH- per mole of Ba(OH)2.
    moles_oh_strong = 2 * moles_baoh2

    # Step 3: Perform the neutralization reaction between the strong acid and strong base.
    # H+ + OH- -> H2O
    # The reactant with fewer moles is the limiting reactant.
    if moles_h_strong > moles_oh_strong:
        moles_h_remaining = moles_h_strong - moles_oh_strong
        moles_oh_remaining = 0
    else:
        moles_oh_remaining = moles_oh_strong - moles_h_strong
        moles_h_remaining = 0

    # Step 4: React the remaining strong species with the weak acid.
    # In this case, we have excess OH- (a strong base) reacting with CH3COOH (a weak acid).
    # CH3COOH + OH- -> CH3COO- + H2O
    # Since OH- is a strong base, this reaction goes to completion.
    if moles_oh_remaining > moles_ch3cooh:
        # The strong base is in excess, so it will determine the final pH.
        final_moles_oh = moles_oh_remaining - moles_ch3cooh
        final_moles_h = 0
    else:
        # This case would result in a buffer or a solution of only conjugate base,
        # but based on the numbers, it's not the scenario here.
        # We can set final_moles_oh to 0 to indicate an error if this path is taken.
        final_moles_oh = 0 # This path should not be taken for this problem.

    # Step 5: Calculate the final pH.
    # The pH is determined by the concentration of the excess strong base.
    total_volume = vol_ch3cooh + vol_hcl + vol_baoh2
    
    if final_moles_oh > 0:
        final_conc_oh = final_moles_oh / total_volume
        poh = -math.log10(final_conc_oh)
        calculated_ph = 14 - poh
    else:
        # This would indicate a different chemical scenario (e.g., buffer, excess acid)
        # or an error in the logic for the given inputs.
        return "Calculation error: The final solution is not dominated by excess strong base as expected."

    # Step 6: Compare the calculated pH with the LLM's answer.
    # A small tolerance is used to account for floating-point and rounding differences.
    tolerance = 0.01
    if abs(calculated_ph - llm_answer_ph) < tolerance:
        return "Correct"
    else:
        reason = (
            f"The provided answer is incorrect.\n"
            f"The calculation shows the following:\n"
            f"1. Moles of H+ from strong acid (HCl): {moles_h_strong:.3f} mol.\n"
            f"2. Moles of OH- from strong base (Ba(OH)2): 2 * (0.3 L * 0.3 M) = {moles_oh_strong:.3f} mol.\n"
            f"3. After strong acid-base neutralization, excess OH- is {moles_oh_strong:.3f} - {moles_h_strong:.3f} = {moles_oh_remaining:.3f} mol.\n"
            f"4. Moles of weak acid (CH3COOH): {moles_ch3cooh:.3f} mol.\n"
            f"5. The excess OH- reacts with the weak acid. Final excess OH- is {moles_oh_remaining:.3f} - {moles_ch3cooh:.3f} = {final_moles_oh:.3f} mol.\n"
            f"6. Total volume is {total_volume:.1f} L.\n"
            f"7. Final [OH-] = {final_moles_oh:.3f} mol / {total_volume:.1f} L = {final_conc_oh:.4f} M.\n"
            f"8. pOH = -log10({final_conc_oh:.4f}) = {-math.log10(final_conc_oh):.2f}.\n"
            f"9. pH = 14 - pOH = {calculated_ph:.2f}.\n"
            f"The calculated pH is {calculated_ph:.2f}, which does not match the provided answer of {llm_answer_ph}."
        )
        return reason

# Execute the check and print the result.
# print(check_correctness())