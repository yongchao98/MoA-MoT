import math

def check_answer():
    """
    This function calculates the concentration of orthophosphate ions based on the problem's data
    and checks if the provided answer is correct.
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
        "B": 2.81e-7,
        "C": 5.48e-7,
        "D": 3.97e-7
    }
    
    # The final answer given by the LLM
    llm_answer_key = "A"

    # --- Step 1: Calculate Moles and Concentrations ---
    # Convert volume to Liters
    volume_L = volume_cm3 / 1000.0

    # Moles of KH2PO4 (source of H2PO4-)
    moles_H2PO4_minus = mass_KH2PO4 / mw_KH2PO4

    # Moles of Na2HPO4●2H2O (source of HPO4^2-)
    moles_HPO4_2minus = mass_Na2HPO4_2H2O / mw_Na2HPO4_2H2O

    # Concentration of H2PO4-
    conc_H2PO4_minus = moles_H2PO4_minus / volume_L

    # Concentration of HPO4^2-
    conc_HPO4_2minus = moles_HPO4_2minus / volume_L

    # --- Step 2: Calculate Hydrogen Ion Concentration [H+] ---
    # The buffer system is H2PO4- / HPO4^2-, governed by Ka2.
    # Ka2 = [H+][HPO4^2-] / [H2PO4-]
    # [H+] = Ka2 * ([H2PO4-] / [HPO4^2-])
    conc_H_plus = Ka2 * (conc_H2PO4_minus / conc_HPO4_2minus)

    # --- Step 3: Calculate Orthophosphate Ion Concentration [PO4^3-] ---
    # The formation of PO4^3- is governed by Ka3.
    # Ka3 = [H+][PO4^3-] / [HPO4^2-]
    # [PO4^3-] = Ka3 * ([HPO4^2-] / [H+])
    calculated_conc_PO4_3minus = Ka3 * (conc_HPO4_2minus / conc_H_plus)

    # --- Step 4: Check the correctness of the LLM's answer ---
    llm_answer_value = options.get(llm_answer_key)
    if llm_answer_value is None:
        return f"The provided answer key '{llm_answer_key}' is not a valid option (A, B, C, D)."

    # Compare the calculated value with the LLM's chosen option value
    # We use math.isclose for robust floating-point comparison.
    # A relative tolerance of 1% (0.01) is reasonable for this kind of chemistry problem
    # due to rounding in the options.
    if math.isclose(calculated_conc_PO4_3minus, llm_answer_value, rel_tol=0.01):
        return "Correct"
    else:
        # Find which option is the correct one
        correct_key = None
        for key, value in options.items():
            if math.isclose(calculated_conc_PO4_3minus, value, rel_tol=0.01):
                correct_key = key
                break
        
        reason = (f"The answer is incorrect.\n"
                  f"The calculated concentration of orthophosphate ions is approximately {calculated_conc_PO4_3minus:.3e} M.\n"
                  f"This value matches option {correct_key} ({options.get(correct_key):.3e} M).\n"
                  f"The provided answer was {llm_answer_key} ({llm_answer_value:.3e} M), which is not correct.")
        return reason

# Run the check
result = check_answer()
print(result)