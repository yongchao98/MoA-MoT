import math

def check_answer():
    """
    This function checks the correctness of the provided answer for the chemistry problem.
    It recalculates the concentration of orthophosphate ions based on the given data and compares it to the answer.
    """
    # Given values from the question
    volume_cm3 = 200.00
    mass_KH2PO4 = 1.00  # g
    mw_KH2PO4 = 136.09  # g/mol
    mass_Na2HPO4_2H2O = 1.00  # g
    mw_Na2HPO4_2H2O = 177.99  # g/mol
    Ka2 = 6.2e-8
    Ka3 = 1.8e-12

    # The options provided in the question
    options = {
        "A": 6.24e-7,
        "B": 5.48e-7,
        "C": 3.97e-7,
        "D": 2.81e-7
    }
    
    # The final answer provided by the LLM
    llm_answer_letter = "A"

    # --- Step-by-step calculation ---

    # 1. Calculate moles of buffer components
    moles_H2PO4_minus = mass_KH2PO4 / mw_KH2PO4
    moles_HPO4_2minus = mass_Na2HPO4_2H2O / mw_Na2HPO4_2H2O

    # 2. Calculate initial concentrations
    volume_L = volume_cm3 / 1000.0
    conc_H2PO4_minus = moles_H2PO4_minus / volume_L
    conc_HPO4_2minus = moles_HPO4_2minus / volume_L

    # 3. Calculate the hydrogen ion concentration [H+] using Ka2
    # Ka2 = [H+][HPO4^2-] / [H2PO4^-]
    # [H+] = Ka2 * [H2PO4^-] / [HPO4^2-]
    conc_H_plus = Ka2 * (conc_H2PO4_minus / conc_HPO4_2minus)

    # 4. Calculate the orthophosphate ion concentration [PO4^3-] using Ka3
    # Ka3 = [H+][PO4^3-] / [HPO4^2-]
    # [PO4^3-] = Ka3 * [HPO4^2-] / [H+]
    calculated_conc_PO4_3minus = Ka3 * (conc_HPO4_2minus / conc_H_plus)

    # --- Verification ---
    
    # Get the value corresponding to the LLM's chosen answer letter
    llm_answer_value = options.get(llm_answer_letter)

    if llm_answer_value is None:
        return f"Invalid answer option '{llm_answer_letter}'. Valid options are A, B, C, D."

    # Check if the calculated value is close to the value of the chosen option.
    # A relative tolerance of 1% is reasonable given the rounding in the options.
    if math.isclose(calculated_conc_PO4_3minus, llm_answer_value, rel_tol=0.01):
        return "Correct"
    else:
        # Find which option the calculation actually matches
        correct_option = None
        for option, value in options.items():
            if math.isclose(calculated_conc_PO4_3minus, value, rel_tol=0.01):
                correct_option = option
                break
        
        reason = (f"The calculated concentration of orthophosphate ions is approximately {calculated_conc_PO4_3minus:.3e} M. "
                  f"The provided answer is '{llm_answer_letter}', which corresponds to a value of {llm_answer_value:.3e} M. "
                  f"The calculation does not support the provided answer. The calculated value matches option '{correct_option}'.")
        return reason

# Run the check
result = check_answer()
print(result)