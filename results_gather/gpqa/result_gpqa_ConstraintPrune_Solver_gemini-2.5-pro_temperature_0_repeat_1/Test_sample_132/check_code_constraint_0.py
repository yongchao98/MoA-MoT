import math

def check_answer():
    """
    This function checks the correctness of the LLM's answer to the chemistry problem.
    """
    # --- Problem Constraints and Given Values ---
    volume_L = 200.00 / 1000.0  # in Liters
    mass_kh2po4_g = 1.00
    mw_kh2po4 = 136.09  # g/mol
    mass_na2hpo4_2h2o_g = 1.00
    mw_na2hpo4_2h2o = 177.99  # g/mol
    ka2 = 6.2e-8
    ka3 = 1.8e-12

    # The LLM's answer is A, which corresponds to this value.
    llm_answer_value = 6.24e-7  # M

    # --- Step 1: Calculate initial molar concentrations of buffer components ---
    # Moles of KH2PO4 provides moles of H2PO4-
    moles_h2po4_minus = mass_kh2po4_g / mw_kh2po4
    conc_h2po4_minus = moles_h2po4_minus / volume_L

    # Moles of Na2HPO4.2H2O provides moles of HPO4^2-
    moles_hpo4_2minus = mass_na2hpo4_2h2o_g / mw_na2hpo4_2h2o
    conc_hpo4_2minus = moles_hpo4_2minus / volume_L

    # --- Step 2: Calculate the hydrogen ion concentration [H+] ---
    # The solution is a buffer, so we use the Ka2 equilibrium expression.
    # H2PO4- <=> H+ + HPO4^2-
    # Ka2 = [H+][HPO4^2-] / [H2PO4-]
    # We assume the initial concentrations are a good approximation of equilibrium concentrations.
    # [H+] = Ka2 * [H2PO4-] / [HPO4^2-]
    conc_h_plus = ka2 * (conc_h2po4_minus / conc_hpo4_2minus)

    # --- Step 3: Calculate the orthophosphate ion concentration [PO4^3-] ---
    # This is governed by the Ka3 equilibrium.
    # HPO4^2- <=> H+ + PO4^3-
    # Ka3 = [H+][PO4^3-] / [HPO4^2-]
    # [PO4^3-] = Ka3 * [HPO4^2-] / [H+]
    calculated_conc_po4_3minus = ka3 * (conc_hpo4_2minus / conc_h_plus)

    # --- Final Verification ---
    # Check if the calculated result is close to the provided answer's value.
    # A relative tolerance of 1% (1e-2) is sufficient given the separation of the options.
    if math.isclose(calculated_conc_po4_3minus, llm_answer_value, rel_tol=1e-2):
        return "Correct"
    else:
        return (
            f"Incorrect.\n"
            f"The provided answer is {llm_answer_value:.2e} M (Option A).\n"
            f"However, the calculated value is {calculated_conc_po4_3minus:.2e} M.\n"
            f"Calculation breakdown:\n"
            f"1. [H₂PO₄⁻] = (1.00 g / 136.09 g/mol) / 0.200 L = {conc_h2po4_minus:.5f} M\n"
            f"2. [HPO₄²⁻] = (1.00 g / 177.99 g/mol) / 0.200 L = {conc_hpo4_2minus:.5f} M\n"
            f"3. [H⁺] = Ka₂ * ([H₂PO₄⁻] / [HPO₄²⁻]) = {conc_h_plus:.2e} M\n"
            f"4. [PO₄³⁻] = Ka₃ * ([HPO₄²⁻] / [H⁺]) = {calculated_conc_po4_3minus:.2e} M\n"
            f"The calculated value {calculated_conc_po4_3minus:.2e} M is significantly different from the answer {llm_answer_value:.2e} M."
        )

# Run the check and print the result
print(check_answer())