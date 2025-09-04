import math

def check_ph_calculation():
    """
    This function checks the correctness of the provided answer for a pH calculation problem.
    It recalculates the pH from the initial conditions and compares it to the given answer.
    """
    # --- Problem Definition ---
    # Solution mixture:
    # 1. 500 mL of 0.1 M CH3COOH (weak acid)
    # 2. 400 mL of 0.2 M HCl (strong acid)
    # 3. 300 mL of 0.3 M Ba(OH)2 (strong base)
    # The provided answer is C, which corresponds to a pH of 12.62.

    # --- Step 1: Calculate initial moles of all species ---
    # Moles = Molarity (mol/L) * Volume (L)
    moles_ch3cooh = 0.1 * 0.500  # 0.05 mol
    moles_hcl = 0.2 * 0.400      # 0.08 mol
    moles_baoh2 = 0.3 * 0.300    # 0.09 mol

    # --- Step 2: Calculate total moles of acidic and basic ions ---
    # HCl is a strong acid, dissociates completely: HCl -> H+ + Cl-
    moles_h_plus = moles_hcl

    # Ba(OH)2 is a strong base, dissociates completely: Ba(OH)2 -> Ba^2+ + 2OH-
    # Note the factor of 2 for OH- ions.
    moles_oh_minus = moles_baoh2 * 2  # 0.09 * 2 = 0.18 mol

    # CH3COOH is a weak acid.
    moles_weak_acid = moles_ch3cooh

    # --- Step 3: Perform neutralization reactions (stoichiometry) ---
    # Priority 1: Strong acid reacts with strong base.
    # H+ + OH- -> H2O
    if moles_h_plus >= moles_oh_minus:
        # This case is not applicable here, but included for completeness
        remaining_h_plus = moles_h_plus - moles_oh_minus
        remaining_oh_minus = 0
    else:
        remaining_oh_minus = moles_oh_minus - moles_h_plus
        remaining_h_plus = 0
    # Calculation: remaining_oh_minus = 0.18 - 0.08 = 0.10 mol

    # Priority 2: Remaining strong species (OH-) reacts with weak acid.
    # CH3COOH + OH- -> CH3COO- + H2O
    if remaining_oh_minus > 0:
        if moles_weak_acid >= remaining_oh_minus:
            # Not applicable here
            final_moles_weak_acid = moles_weak_acid - remaining_oh_minus
            final_moles_oh = 0
        else:
            final_moles_oh = remaining_oh_minus - moles_weak_acid
            final_moles_weak_acid = 0
        # Calculation: final_moles_oh = 0.10 - 0.05 = 0.05 mol
    else:
        # This would mean excess acid, not applicable here.
        final_moles_oh = 0

    # --- Step 4: Calculate final pH ---
    # The final pH is determined by the concentration of the excess strong base (OH-).
    # The contribution from the hydrolysis of the conjugate base (CH3COO-) is negligible.
    total_volume = 0.500 + 0.400 + 0.300  # 1.2 L

    if final_moles_oh <= 0:
        return "Calculation Error: The final solution should have excess OH-, but the calculation resulted in zero or negative moles."

    # Calculate the concentration of OH-
    concentration_oh = final_moles_oh / total_volume
    # Calculation: concentration_oh = 0.05 mol / 1.2 L = 0.041667 M

    # Calculate pOH and then pH
    poh = -math.log10(concentration_oh)
    calculated_ph = 14.0 - poh
    # Calculation: pOH = 1.38, pH = 14 - 1.38 = 12.62

    # --- Step 5: Verify the answer ---
    expected_ph = 12.62  # From option C
    # Check if the calculated pH matches the expected pH within a small tolerance
    if abs(calculated_ph - expected_ph) < 0.01:
        return "Correct"
    else:
        return (f"Incorrect. The provided answer is {expected_ph}, but the calculated pH is {calculated_ph:.2f}. "
                f"The final moles of excess OH- were calculated to be {final_moles_oh:.4f} in a total volume of {total_volume} L, "
                f"leading to a different pH.")

# Execute the check and print the result
result = check_ph_calculation()
print(result)