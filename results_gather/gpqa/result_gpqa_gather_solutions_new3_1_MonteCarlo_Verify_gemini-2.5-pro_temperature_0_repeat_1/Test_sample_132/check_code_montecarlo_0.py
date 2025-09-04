import math

def check_answer_correctness():
    """
    This function checks the correctness of the provided answer to the chemistry question.
    It recalculates the concentration of orthophosphate ions based on the problem statement
    and compares it to the value of the chosen option.
    """
    
    # --- Define constants and given values from the question ---
    volume_cm3 = 200.00
    mass_kh2po4 = 1.00  # g
    mw_kh2po4 = 136.09  # g/mol
    mass_na2hpo4_2h2o = 1.00  # g
    mw_na2hpo4_2h2o = 177.99  # g/mol
    ka2 = 6.2e-8
    ka3 = 1.8e-12

    # --- Define the options from the question ---
    options = {
        "A": 3.97e-7,
        "B": 5.48e-7,
        "C": 2.81e-7,
        "D": 6.24e-7
    }
    
    # --- The final answer provided by the LLM to be checked ---
    # The LLM's final answer is <<<D>>>
    llm_answer_key = "D"
    
    # --- Step 1: Calculate moles of each component ---
    # Moles of KH2PO4 provides moles of H2PO4-
    moles_h2po4_minus = mass_kh2po4 / mw_kh2po4
    # Moles of Na2HPO4●2H2O provides moles of HPO4^2-
    moles_hpo4_2minus = mass_na2hpo4_2h2o / mw_na2hpo4_2h2o

    # --- Step 2: Calculate initial concentrations ---
    volume_L = volume_cm3 / 1000.0
    conc_h2po4_minus = moles_h2po4_minus / volume_L
    conc_hpo4_2minus = moles_hpo4_2minus / volume_L

    # --- Step 3: Calculate the hydrogen ion concentration [H+] ---
    # The buffer system is H2PO4- / HPO4^2-, governed by Ka2.
    # Equilibrium: H2PO4- <=> H+ + HPO4^2-
    # Ka2 = [H+][HPO4^2-] / [H2PO4-]
    # Rearranging gives: [H+] = Ka2 * [H2PO4-] / [HPO4^2-]
    conc_h_plus = ka2 * (conc_h2po4_minus / conc_hpo4_2minus)

    # --- Step 4: Calculate the orthophosphate ion concentration [PO4^3-] ---
    # The formation of PO4^3- is governed by Ka3.
    # Equilibrium: HPO4^2- <=> H+ + PO4^3-
    # Ka3 = [H+][PO4^3-] / [HPO4^2-]
    # Rearranging gives: [PO4^3-] = Ka3 * [HPO4^2-] / [H+]
    calculated_conc_po4_3minus = ka3 * (conc_hpo4_2minus / conc_h_plus)

    # --- Step 5: Check the correctness of the LLM's answer ---
    if llm_answer_key not in options:
        return f"Invalid answer key '{llm_answer_key}'. The key must be one of {list(options.keys())}."

    expected_value = options[llm_answer_key]

    # Compare the calculated value with the value from the chosen option.
    # A relative tolerance of 1% is used to account for potential rounding differences in the options.
    if math.isclose(calculated_conc_po4_3minus, expected_value, rel_tol=0.01):
        return "Correct"
    else:
        # Find the option that is numerically closest to the calculated answer
        closest_option_key = min(options, key=lambda k: abs(options[k] - calculated_conc_po4_3minus))
        
        reason = (
            f"The provided answer '{llm_answer_key}' is incorrect.\n"
            f"The value corresponding to option {llm_answer_key} is {expected_value:.2e} M.\n"
            f"The correct calculation is as follows:\n"
            f"1. [H₂PO₄⁻] = (1.00 g / 136.09 g/mol) / 0.200 L = {conc_h2po4_minus:.4f} M\n"
            f"2. [HPO₄²⁻] = (1.00 g / 177.99 g/mol) / 0.200 L = {conc_hpo4_2minus:.4f} M\n"
            f"3. [H⁺] = Ka₂ * ([H₂PO₄⁻]/[HPO₄²⁻]) = {conc_h_plus:.3e} M\n"
            f"4. [PO₄³⁻] = Ka₃ * ([HPO₄²⁻]/[H⁺]) = {calculated_conc_po4_3minus:.3e} M\n"
            f"The calculated result {calculated_conc_po4_3minus:.3e} M is closest to option '{closest_option_key}' ({options[closest_option_key]:.2e} M), not option '{llm_answer_key}'."
        )
        return reason

# Execute the check and print the result
print(check_answer_correctness())