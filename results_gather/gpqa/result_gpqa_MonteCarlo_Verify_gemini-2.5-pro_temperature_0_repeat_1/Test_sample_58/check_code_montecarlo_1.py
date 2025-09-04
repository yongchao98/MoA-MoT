import math

def check_ph_correctness():
    """
    This function checks the correctness of the given answer for a pH calculation problem.
    It calculates the pH of a solution formed by mixing a weak acid, a strong acid, and a strong base.
    
    The problem involves:
    - 500 mL of 0.1 M CH3COOH (weak acid)
    - 400 mL of 0.2 M HCl (strong acid)
    - 300 mL of 0.3 M Ba(OH)2 (strong base)
    
    The function compares its calculated pH to the provided answer (12.62) and returns
    "Correct" or a detailed explanation of the error.
    """
    
    # --- Given Data from the Question ---
    # Component 1: Acetic Acid (CH3COOH), a weak acid (HA)
    vol_ha = 0.500  # L (500 mL)
    conc_ha = 0.1   # M

    # Component 2: Hydrochloric Acid (HCl), a strong acid
    vol_hcl = 0.400 # L (400 mL)
    conc_hcl = 0.2  # M

    # Component 3: Barium Hydroxide (Ba(OH)2), a strong base
    vol_baoh2 = 0.300 # L (300 mL)
    conc_baoh2 = 0.3  # M

    # The answer from the other LLM to be checked
    given_answer_ph = 12.62

    # --- Step 1: Calculate initial moles of all reactive species ---
    # Moles = Molarity * Volume
    
    # Moles of weak acid
    moles_ha = conc_ha * vol_ha
    
    # Moles of H+ from strong acid (HCl -> H+ + Cl-)
    moles_h_strong = conc_hcl * vol_hcl
    
    # Moles of OH- from strong base (Ba(OH)2 -> Ba^2+ + 2OH-)
    # A critical constraint is that 1 mole of Ba(OH)2 produces 2 moles of OH-
    moles_oh_strong = conc_baoh2 * vol_baoh2 * 2

    # --- Step 2: Perform neutralization reactions in order of strength ---
    # First, the strong acid reacts completely with the strong base.
    # H+ + OH- -> H2O
    moles_oh_after_strong_reaction = moles_oh_strong - moles_h_strong

    # Second, the remaining strong species (in this case, OH-) reacts with the weak acid.
    # CH3COOH + OH- -> CH3COO- + H2O
    final_moles_oh = moles_oh_after_strong_reaction - moles_ha

    # --- Step 3: Calculate final pH ---
    # The final pH is determined by the excess species. In this case, there is an
    # excess of strong base (OH-), which will dominate the pH of the solution.
    
    # Calculate total volume of the solution
    total_volume = vol_ha + vol_hcl + vol_baoh2

    # Check if there is indeed an excess of strong base
    if final_moles_oh > 0:
        # Calculate the final concentration of OH-
        final_conc_oh = final_moles_oh / total_volume
        
        # Calculate pOH from [OH-]
        poh = -math.log10(final_conc_oh)
        
        # Calculate pH from pOH (pH + pOH = 14)
        calculated_ph = 14.0 - poh
    else:
        # This block would handle cases with excess acid or a buffer solution,
        # but based on the initial numbers, we expect excess strong base.
        # If the calculation leads here, it indicates a different chemical scenario.
        return "The calculation did not result in an excess of strong base as expected."

    # --- Step 4: Compare calculated pH with the given answer ---
    # We use math.isclose() to account for potential floating-point inaccuracies.
    # A relative tolerance of 1% is suitable for this kind of chemical calculation.
    if math.isclose(calculated_ph, given_answer_ph, rel_tol=1e-2):
        return "Correct"
    else:
        # If the answer is incorrect, provide a detailed breakdown of the correct calculation.
        reason = (
            f"The provided answer {given_answer_ph} is incorrect. The calculated pH is {calculated_ph:.2f}.\n\n"
            f"Here is the step-by-step calculation:\n"
            f"1. Calculate Initial Moles:\n"
            f"   - Moles CH3COOH = {conc_ha} M * {vol_ha} L = {moles_ha:.4f} mol\n"
            f"   - Moles H+ (from HCl) = {conc_hcl} M * {vol_hcl} L = {moles_h_strong:.4f} mol\n"
            f"   - Moles OH- (from Ba(OH)2) = {conc_baoh2} M * {vol_baoh2} L * 2 = {moles_oh_strong:.4f} mol. (Constraint: Ba(OH)2 provides 2 OH- ions)\n\n"
            f"2. Neutralize Strong Acid with Strong Base:\n"
            f"   - Total moles of base ({moles_oh_strong:.4f}) is greater than moles of strong acid ({moles_h_strong:.4f}).\n"
            f"   - Moles of OH- remaining = {moles_oh_strong:.4f} - {moles_h_strong:.4f} = {moles_oh_after_strong_reaction:.4f} mol\n\n"
            f"3. Neutralize Weak Acid with Remaining Strong Base:\n"
            f"   - Remaining moles of OH- ({moles_oh_after_strong_reaction:.4f}) is greater than moles of weak acid ({moles_ha:.4f}).\n"
            f"   - Final moles of excess OH- = {moles_oh_after_strong_reaction:.4f} - {moles_ha:.4f} = {final_moles_oh:.4f} mol\n\n"
            f"4. Calculate Final pH:\n"
            f"   - Total Volume = {vol_ha} L + {vol_hcl} L + {vol_baoh2} L = {total_volume:.1f} L\n"
            f"   - Final [OH-] = {final_moles_oh:.4f} mol / {total_volume:.1f} L = {final_conc_oh:.5f} M\n"
            f"   - pOH = -log10({final_conc_oh:.5f}) = {-math.log10(final_conc_oh):.2f}\n"
            f"   - pH = 14 - pOH = 14 - {-math.log10(final_conc_oh):.2f} = {calculated_ph:.2f}\n"
        )
        return reason

# Execute the checking function and print its output.
print(check_ph_correctness())