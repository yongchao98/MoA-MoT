import math

def check_ph_calculation():
    """
    This function checks the correctness of the pH calculation for the given mixture.
    It recalculates the pH based on the problem's constraints and compares it to the provided answer.
    """
    # --- Define initial conditions from the question ---
    # Weak acid: 500 mL of 0.1 M CH3COOH
    vol_ch3cooh = 0.500  # L
    conc_ch3cooh = 0.1   # M

    # Strong acid: 400 mL of 0.2 M HCl
    vol_hcl = 0.400     # L
    conc_hcl = 0.2       # M

    # Strong base: 300 mL of 0.3 M Ba(OH)2
    vol_baoh2 = 0.300   # L
    conc_baoh2 = 0.3     # M

    # The answer to check
    expected_ph = 12.62

    # --- Step 1: Calculate the initial moles of each species ---
    # Moles = Concentration * Volume
    moles_ch3cooh = conc_ch3cooh * vol_ch3cooh
    moles_hcl = conc_hcl * vol_hcl
    moles_baoh2 = conc_baoh2 * vol_baoh2

    # --- Step 2: Calculate the total moles of acidic protons (H+) and hydroxide ions (OH-) ---
    # HCl is a strong acid, so it fully dissociates.
    moles_h_from_hcl = moles_hcl
    # CH3COOH is a weak acid, but for neutralization with a strong base, we assume it reacts completely.
    moles_h_from_ch3cooh = moles_ch3cooh
    # Total moles of acid species available for reaction.
    total_moles_acid = moles_h_from_hcl + moles_h_from_ch3cooh

    # Ba(OH)2 is a strong base that provides two OH- ions per formula unit.
    moles_oh_from_baoh2 = 2 * moles_baoh2

    # --- Step 3: Perform the neutralization reaction to find the excess species ---
    # The reaction is H+ + OH- -> H2O. We compare the total moles of acid and base.
    if moles_oh_from_baoh2 > total_moles_acid:
        # Base is in excess
        excess_moles = moles_oh_from_baoh2 - total_moles_acid
        excess_species = "OH-"
    elif total_moles_acid > moles_oh_from_baoh2:
        # Acid is in excess
        excess_moles = total_moles_acid - moles_oh_from_baoh2
        excess_species = "H+"
    else:
        # Perfect neutralization. The pH would be determined by the salt (CH3COOBa), which is basic.
        # However, the problem values make this outcome highly unlikely.
        # For simplicity, we'll assume pH 7 if this happens.
        calculated_ph = 7.0
        if not math.isclose(calculated_ph, expected_ph, abs_tol=0.01):
            return f"Calculation results in perfect neutralization (pH=7.0), which contradicts the expected pH of {expected_ph}."
        else:
            return "Correct"

    # --- Step 4: Calculate the total volume of the final solution ---
    total_volume = vol_ch3cooh + vol_hcl + vol_baoh2

    # --- Step 5: Calculate the final concentration of the excess species ---
    final_concentration = excess_moles / total_volume

    # --- Step 6: Calculate the final pH ---
    if excess_species == "OH-":
        # The presence of a strong base (excess OH-) suppresses the hydrolysis of the weak conjugate base (CH3COO-).
        # Therefore, the pH is determined solely by the concentration of the excess strong base.
        poh = -math.log10(final_concentration)
        calculated_ph = 14.0 - poh
    elif excess_species == "H+":
        # If H+ were in excess, the solution would contain a strong acid (excess H+) and a weak acid (CH3COOH).
        # The pH would be dominated by the concentration of the strong acid.
        calculated_ph = -math.log10(final_concentration)

    # --- Step 7: Compare the calculated pH with the expected answer ---
    # We use a tolerance (e.g., 0.01) for comparing floating-point numbers.
    if math.isclose(calculated_ph, expected_ph, abs_tol=0.01):
        return "Correct"
    else:
        # If the answer is incorrect, provide a step-by-step reason.
        reason = (
            f"The calculated pH of {calculated_ph:.2f} does not match the expected pH of {expected_ph}.\n"
            f"Here is the step-by-step calculation:\n"
            f"1. Moles of total acid species = (moles H+ from HCl) + (moles CH3COOH)\n"
            f"   = ({conc_hcl} M * {vol_hcl} L) + ({conc_ch3cooh} M * {vol_ch3cooh} L)\n"
            f"   = {moles_hcl:.3f} + {moles_ch3cooh:.3f} = {total_moles_acid:.3f} mol\n"
            f"2. Moles of OH- = 2 * (moles Ba(OH)2)\n"
            f"   = 2 * ({conc_baoh2} M * {vol_baoh2} L) = {moles_oh_from_baoh2:.3f} mol\n"
            f"3. Since moles OH- ({moles_oh_from_baoh2:.3f}) > moles acid ({total_moles_acid:.3f}), the base is in excess.\n"
            f"   Moles of excess OH- = {moles_oh_from_baoh2:.3f} - {total_moles_acid:.3f} = {excess_moles:.3f} mol\n"
            f"4. Total volume = {vol_ch3cooh} L + {vol_hcl} L + {vol_baoh2} L = {total_volume:.1f} L\n"
            f"5. Final [OH-] = {excess_moles:.3f} mol / {total_volume:.1f} L = {final_concentration:.5f} M\n"
            f"6. pOH = -log10({final_concentration:.5f}) = {-math.log10(final_concentration):.2f}\n"
            f"7. Calculated pH = 14 - pOH = 14 - {-math.log10(final_concentration):.2f} = {calculated_ph:.2f}\n"
            f"The final calculated pH is {calculated_ph:.2f}, which does not match the provided answer."
        )
        return reason

# Execute the check
result = check_ph_calculation()
print(result)