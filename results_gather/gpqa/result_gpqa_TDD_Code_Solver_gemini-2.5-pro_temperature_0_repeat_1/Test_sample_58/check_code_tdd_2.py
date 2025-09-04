import math

def check_ph_calculation():
    """
    This function checks the correctness of a given answer for a pH calculation problem.
    It calculates the pH of a solution made by mixing a weak acid, a strong acid, and a strong base.
    """
    # --- Problem Data ---
    # Acetic Acid (CH3COOH) - Weak Acid
    vol_ch3cooh = 0.500  # L (500 mL)
    conc_ch3cooh = 0.1   # M

    # Hydrochloric Acid (HCl) - Strong Acid
    vol_hcl = 0.400  # L (400 mL)
    conc_hcl = 0.2   # M

    # Barium Hydroxide (Ba(OH)2) - Strong Base
    vol_baoh2 = 0.300  # L (300 mL)
    conc_baoh2 = 0.3   # M

    # The answer to be checked, corresponding to option C
    given_answer_ph = 12.62

    # --- Step 1: Calculate initial moles of reactive species ---
    # Moles = Concentration * Volume
    
    # Moles of weak acid
    moles_ch3cooh = conc_ch3cooh * vol_ch3cooh
    
    # Moles of H+ from strong acid (HCl -> H+ + Cl-)
    moles_h_from_hcl = conc_hcl * vol_hcl
    
    # Moles of OH- from strong base (Ba(OH)2 -> Ba^2+ + 2OH-)
    # Note the stoichiometry: 1 mole of Ba(OH)2 produces 2 moles of OH-
    moles_oh_from_baoh2 = (conc_baoh2 * vol_baoh2) * 2

    # --- Step 2: Neutralization of strong acid and strong base ---
    # The H+ from the strong acid reacts completely with the OH- from the strong base first.
    # H+ + OH- -> H2O
    
    moles_oh_after_strong_neut = moles_oh_from_baoh2 - moles_h_from_hcl
    
    # --- Step 3: Neutralization involving the weak acid ---
    # The remaining strong base (OH-) reacts with the weak acid (CH3COOH).
    # CH3COOH + OH- -> CH3COO- + H2O
    
    # The limiting reactant is the one with fewer moles.
    # In this case, moles_oh_after_strong_neut (0.10 mol) > moles_ch3cooh (0.05 mol),
    # so all the weak acid is consumed.
    
    final_moles_oh_excess = moles_oh_after_strong_neut - moles_ch3cooh

    # --- Step 4: Calculate final concentration and pH ---
    # The final pH is determined by the concentration of the excess strong base (OH-).
    # The contribution from the hydrolysis of the conjugate base (CH3COO-) is negligible
    # in the presence of a much stronger base.

    # Calculate the total volume of the final solution
    total_volume = vol_ch3cooh + vol_hcl + vol_baoh2

    # Calculate the final concentration of OH-
    final_conc_oh = final_moles_oh_excess / total_volume

    # Calculate pOH from [OH-]
    try:
        poh = -math.log10(final_conc_oh)
    except ValueError:
        return "Calculation Error: Cannot take the log of a non-positive concentration."

    # Calculate pH from pOH (pH + pOH = 14 at 25°C)
    calculated_ph = 14.0 - poh

    # --- Step 5: Compare calculated pH with the given answer ---
    # Use a small tolerance for floating-point comparison.
    if math.isclose(calculated_ph, given_answer_ph, rel_tol=1e-2, abs_tol=1e-2):
        return "Correct"
    else:
        reason = (
            f"The provided answer's pH is {given_answer_ph}, but the calculated pH is {calculated_ph:.2f}.\n"
            f"Constraint check failed.\n"
            f"--- Calculation Breakdown ---\n"
            f"1. Initial Moles:\n"
            f"   - Moles H+ (from strong acid HCl) = {moles_h_from_hcl:.3f} mol\n"
            f"   - Moles OH- (from strong base Ba(OH)2) = {moles_oh_from_baoh2:.3f} mol\n"
            f"   - Moles CH3COOH (weak acid) = {moles_ch3cooh:.3f} mol\n"
            f"2. Neutralization:\n"
            f"   - After strong-strong reaction, moles of OH- remaining = {moles_oh_from_baoh2:.3f} - {moles_h_from_hcl:.3f} = {moles_oh_after_strong_neut:.3f} mol\n"
            f"   - After reacting with weak acid, final excess moles of OH- = {moles_oh_after_strong_neut:.3f} - {moles_ch3cooh:.3f} = {final_moles_oh_excess:.3f} mol\n"
            f"3. Final pH Calculation:\n"
            f"   - Total Volume = {total_volume:.3f} L\n"
            f"   - Final [OH-] = {final_moles_oh_excess:.3f} mol / {total_volume:.3f} L = {final_conc_oh:.4f} M\n"
            f"   - pOH = -log10({final_conc_oh:.4f}) = {poh:.2f}\n"
            f"   - pH = 14 - pOH = {calculated_ph:.2f}\n"
            f"The calculated pH of {calculated_ph:.2f} does not match the given answer of {given_answer_ph}."
        )
        return reason

# Run the check
result = check_ph_calculation()
print(result)