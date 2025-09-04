import math

def check_answer():
    """
    This function checks the correctness of the provided chemistry calculation.
    It recalculates the pH of the mixed solution from scratch and compares it
    to the result given in the answer.
    """
    try:
        # --- Define initial conditions from the question ---
        # Acetic acid (CH3COOH) - Weak Acid
        vol_ch3cooh = 0.500  # L
        conc_ch3cooh = 0.1   # M

        # Hydrochloric acid (HCl) - Strong Acid
        vol_hcl = 0.400      # L
        conc_hcl = 0.2       # M

        # Barium hydroxide (Ba(OH)2) - Strong Base
        vol_baoh2 = 0.300    # L
        conc_baoh2 = 0.3     # M

        # The pH value calculated in the provided answer
        # The answer calculates pH = 12.62, which corresponds to option C.
        expected_ph = 12.62

        # --- Step 1: Calculate initial moles of all reactive species ---
        # Moles of weak acid
        moles_ch3cooh = vol_ch3cooh * conc_ch3cooh

        # Moles of H+ from strong acid
        moles_h_strong = vol_hcl * conc_hcl

        # Moles of OH- from strong base. Note: Ba(OH)2 dissociates into Ba2+ and 2 OH-
        moles_oh_strong = vol_baoh2 * conc_baoh2 * 2

        # --- Step 2: Simulate neutralization reactions ---
        # Reaction 1: Strong acid (H+) reacts with strong base (OH-) first.
        # H+ is the limiting reactant in this step.
        moles_oh_after_hcl = moles_oh_strong - moles_h_strong
        
        # Reaction 2: Remaining strong base (OH-) reacts with the weak acid (CH3COOH).
        # CH3COOH is the limiting reactant in this step.
        moles_oh_final = moles_oh_after_hcl - moles_ch3cooh

        # --- Step 3: Calculate the final pH ---
        # The final solution contains excess strong base (OH-), so its concentration
        # will determine the pH. The contribution from the conjugate base (CH3COO-) is negligible.

        # Calculate the total volume of the solution
        total_volume = vol_ch3cooh + vol_hcl + vol_baoh2

        # Check if there is indeed excess strong base
        if moles_oh_final <= 0:
            return (f"Incorrect calculation logic. The provided answer assumes an excess of strong base (OH-), "
                    f"but my calculation shows this is not the case. Final moles of OH- calculated: {moles_oh_final:.4f}.")

        # Calculate the final concentration of OH-
        final_conc_oh = moles_oh_final / total_volume

        # Calculate pOH and then pH
        pOH = -math.log10(final_conc_oh)
        calculated_ph = 14.0 - pOH

        # --- Step 4: Compare calculated pH with the expected pH ---
        # Use math.isclose() for safe floating-point comparison. A tolerance of 0.01 is reasonable for pH values.
        if math.isclose(calculated_ph, expected_ph, abs_tol=0.01):
            return "Correct"
        else:
            # If the calculation doesn't match, provide the detailed correct calculation.
            reason = (
                f"The answer is incorrect. The final pH is not {expected_ph}.\n"
                f"Here is the correct calculation:\n"
                f"1. Initial Moles:\n"
                f"   - Moles CH3COOH = {moles_ch3cooh:.4f} mol\n"
                f"   - Moles H+ (from HCl) = {moles_h_strong:.4f} mol\n"
                f"   - Moles OH- (from Ba(OH)2) = {moles_oh_strong:.4f} mol\n"
                f"2. Neutralization (Strong Acid vs. Strong Base):\n"
                f"   - Moles OH- remaining = {moles_oh_strong:.4f} - {moles_h_strong:.4f} = {moles_oh_after_hcl:.4f} mol\n"
                f"3. Neutralization (Remaining Base vs. Weak Acid):\n"
                f"   - Final moles of excess OH- = {moles_oh_after_hcl:.4f} - {moles_ch3cooh:.4f} = {moles_oh_final:.4f} mol\n"
                f"4. Final pH Calculation:\n"
                f"   - Total Volume = {total_volume:.3f} L\n"
                f"   - Final [OH-] = {moles_oh_final:.4f} mol / {total_volume:.3f} L = {final_conc_oh:.4f} M\n"
                f"   - pOH = -log10({final_conc_oh:.4f}) = {pOH:.2f}\n"
                f"   - Calculated pH = 14 - {pOH:.2f} = {calculated_ph:.2f}\n"
                f"The calculated pH is {calculated_ph:.2f}, which does not match the expected value of {expected_ph}."
            )
            return reason

    except Exception as e:
        return f"An error occurred during the check: {e}"

# Execute the check and print the result
print(check_answer())