import math

def check_ph_calculation():
    """
    This function checks the correctness of the pH calculation for the given mixture.
    It recalculates the pH based on first principles and compares it to the provided answer.
    """
    # --- Problem Parameters ---
    # Acetic Acid (CH3COOH) - Weak Acid
    vol_acid_weak = 500 / 1000  # L
    conc_acid_weak = 0.1        # M

    # Hydrochloric Acid (HCl) - Strong Acid
    vol_acid_strong = 400 / 1000  # L
    conc_acid_strong = 0.2        # M

    # Barium Hydroxide (Ba(OH)2) - Strong Base
    vol_base_strong = 300 / 1000  # L
    conc_base_strong = 0.3        # M

    # The pH value from the provided answer
    answer_ph = 12.62

    # --- Step 1: Calculate initial moles of each substance ---
    moles_acid_weak = vol_acid_weak * conc_acid_weak
    moles_acid_strong = vol_acid_strong * conc_acid_strong
    moles_base_strong = vol_base_strong * conc_base_strong

    # --- Step 2: Calculate total moles of acidic protons (H+) and hydroxide ions (OH-) ---
    # Total moles of H+ available for reaction from both strong and weak acids.
    # The strong base will react with all of it.
    total_moles_h_plus = moles_acid_strong + moles_acid_weak

    # Total moles of OH- from the strong base.
    # Ba(OH)2 is a dihydroxic base, so it provides 2 moles of OH- per mole of Ba(OH)2.
    total_moles_oh_minus = moles_base_strong * 2

    # --- Step 3: Determine the excess reactant after neutralization ---
    if not math.isclose(total_moles_oh_minus, 0.18):
        return f"Calculation of total moles of OH- is incorrect. Calculated {total_moles_oh_minus}, but should be 0.090 mol * 2 = 0.180 mol."
    
    if not math.isclose(total_moles_h_plus, 0.13):
        return f"Calculation of total moles of H+ is incorrect. Calculated {total_moles_h_plus}, but should be 0.080 mol (from HCl) + 0.050 mol (from CH3COOH) = 0.130 mol."

    if total_moles_oh_minus > total_moles_h_plus:
        # Base is in excess
        excess_moles_oh = total_moles_oh_minus - total_moles_h_plus
    elif total_moles_h_plus > total_moles_oh_minus:
        # This case is not expected based on the numbers, but included for completeness.
        excess_moles_h = total_moles_h_plus - total_moles_oh_minus
        total_volume = vol_acid_weak + vol_acid_strong + vol_base_strong
        final_h_conc = excess_moles_h / total_volume
        calculated_ph = -math.log10(final_h_conc)
        return f"The calculation shows excess acid, leading to a pH of {calculated_ph:.2f}. The provided answer suggests an excess of base."
    else:
        # Neutralization is exact
        return "The calculation shows exact neutralization (pH=7), which contradicts the provided answer."

    # --- Step 4: Calculate the concentration of the excess OH- ---
    total_volume = vol_acid_weak + vol_acid_strong + vol_base_strong
    if not math.isclose(total_volume, 1.2):
        return f"Total volume calculation is incorrect. Calculated {total_volume} L, but should be 0.5 + 0.4 + 0.3 = 1.2 L."

    final_oh_conc = excess_moles_oh / total_volume

    # --- Step 5: Calculate pOH and then pH ---
    # The pH is determined by the excess strong base. The contribution from the
    # conjugate base (CH3COO-) is negligible in the presence of a strong base.
    poh = -math.log10(final_oh_conc)
    calculated_ph = 14.0 - poh

    # --- Step 6: Compare the calculated pH with the provided answer ---
    # We use math.isclose to account for potential floating-point and rounding differences.
    # A tolerance of 0.01 is reasonable for pH values.
    if math.isclose(calculated_ph, answer_ph, abs_tol=0.01):
        return "Correct"
    else:
        return (f"The final pH is incorrect. "
                f"The provided answer is {answer_ph}, but the calculated pH is {calculated_ph:.2f}. "
                f"Calculation details: "
                f"Total H+ moles = {total_moles_h_plus:.3f}. "
                f"Total OH- moles = {total_moles_oh_minus:.3f}. "
                f"Excess OH- moles = {excess_moles_oh:.3f}. "
                f"Total Volume = {total_volume:.3f} L. "
                f"Final [OH-] = {final_oh_conc:.5f} M. "
                f"pOH = -log({final_oh_conc:.5f}) = {poh:.2f}. "
                f"pH = 14 - {poh:.2f} = {calculated_ph:.2f}.")

# Run the check
result = check_ph_calculation()
print(result)