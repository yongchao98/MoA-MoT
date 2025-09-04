import math

def check_ph_calculation():
    """
    This function checks the correctness of the pH calculation for the given chemical mixture.
    It recalculates the pH from scratch and compares it to the provided answer.
    """
    
    # --- Define problem constraints and given values ---
    # 500 mL of 0.1 M CH3COOH
    vol_ch3cooh = 0.500  # L
    conc_ch3cooh = 0.1   # M
    # 400 mL of 0.2 M HCl
    vol_hcl = 0.400      # L
    conc_hcl = 0.2       # M
    # 300 mL of 0.3 M Ba(OH)2
    vol_baoh2 = 0.300    # L
    conc_baoh2 = 0.3     # M

    # The final answer from the LLM to be checked
    llm_answer_value = 12.62
    llm_answer_option = 'C'
    options = {'A': 8.68, 'B': 8.92, 'C': 12.62, 'D': 1.38}

    # --- Step-by-step calculation ---

    # 1. Calculate initial moles of acidic and basic species
    moles_ch3cooh = vol_ch3cooh * conc_ch3cooh
    moles_h_from_hcl = vol_hcl * conc_hcl
    # Ba(OH)2 provides 2 OH- ions per formula unit
    moles_oh_from_baoh2 = vol_baoh2 * conc_baoh2 * 2

    # 2. Determine total moles of acid and base for neutralization
    # The strong base will neutralize all available acidic protons
    total_moles_acid = moles_h_from_hcl + moles_ch3cooh
    total_moles_base = moles_oh_from_baoh2

    # 3. Find the excess reactant
    if total_moles_base > total_moles_acid:
        excess_moles_oh = total_moles_base - total_moles_acid
        
        # 4. Calculate final concentration of excess OH-
        total_volume = vol_ch3cooh + vol_hcl + vol_baoh2
        final_conc_oh = excess_moles_oh / total_volume
        
        # 5. Calculate pOH and then pH
        poh = -math.log10(final_conc_oh)
        calculated_ph = 14.0 - poh
    else:
        # This branch handles cases with excess acid or exact neutralization,
        # which are not the case for this specific problem.
        return "Calculation logic error: The problem was expected to have excess base."

    # --- Verification ---

    # Check 1: Does the calculated pH match the value in the LLM's reasoning?
    if not math.isclose(calculated_ph, llm_answer_value, rel_tol=1e-2):
        return (f"Incorrect calculation. The calculated pH is {calculated_ph:.2f}, "
                f"but the answer states it is {llm_answer_value}.")

    # Check 2: Does the LLM's final option letter match the calculated value?
    if not math.isclose(options.get(llm_answer_option, None), llm_answer_value, rel_tol=1e-2):
        return (f"Incorrect final answer format. The answer choice '{llm_answer_option}' corresponds to pH {options.get(llm_answer_option)}, "
                f"which does not match the calculated value of {llm_answer_value:.2f}.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Execute the check and print the result
result = check_ph_calculation()
print(result)