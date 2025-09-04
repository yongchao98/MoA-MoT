import math

def check_answer():
    """
    This function checks the correctness of the pH calculation for the given mixture.
    """
    # --- 1. Define initial conditions from the question ---
    # CH3COOH (weak acid)
    vol_ch3cooh = 0.500  # L
    conc_ch3cooh = 0.1   # M

    # HCl (strong acid)
    vol_hcl = 0.400      # L
    conc_hcl = 0.2       # M

    # Ba(OH)2 (strong base)
    vol_baoh2 = 0.300    # L
    conc_baoh2 = 0.3     # M

    # --- 2. Calculate initial moles of all species ---
    moles_ch3cooh = vol_ch3cooh * conc_ch3cooh
    moles_h_from_hcl = vol_hcl * conc_hcl
    
    # Constraint: Ba(OH)2 is a dihydroxide base, providing 2 moles of OH- per mole.
    moles_oh_from_baoh2 = vol_baoh2 * conc_baoh2 * 2

    # --- 3. Determine total moles of acid and base and find the excess ---
    total_moles_acid = moles_ch3cooh + moles_h_from_hcl
    total_moles_base = moles_oh_from_baoh2

    # The final solution will be basic as moles of base > moles of acid
    if total_moles_base <= total_moles_acid:
        return "Incorrect: The calculation should result in an excess of base, but it doesn't."

    excess_moles_oh = total_moles_base - total_moles_acid

    # --- 4. Calculate the total volume ---
    total_volume = vol_ch3cooh + vol_hcl + vol_baoh2

    # --- 5. Calculate the final concentration of the excess species ---
    final_conc_oh = excess_moles_oh / total_volume

    # --- 6. Calculate pOH and then pH ---
    # Constraint: For a basic solution, calculate pOH first, then find pH.
    poh = -math.log10(final_conc_oh)
    calculated_ph = 14 - poh

    # --- 7. Compare the calculated pH with the value from the answer (A: 12.62) ---
    answer_ph = 12.62
    
    # Check if the calculated pH is close to the answer's pH value.
    # A small tolerance is used for floating-point comparisons.
    if math.isclose(calculated_ph, answer_ph, rel_tol=1e-3):
        return "Correct"
    else:
        return (f"Incorrect: The calculated pH is {calculated_ph:.2f}, "
                f"but the answer corresponds to a pH of {answer_ph}. "
                f"The provided answer is not consistent with the calculation.")

# Run the check
result = check_answer()
print(result)