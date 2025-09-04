import math

def check_chemistry_answer():
    """
    This function checks the correctness of the provided answer for the chemistry problem.
    It recalculates the concentration of orthophosphate ions based on the given data and compares it
    to the selected option.
    """
    # --- Problem Constraints and Given Data ---
    volume_cm3 = 200.00
    mass_KH2PO4 = 1.00  # g
    mw_KH2PO4 = 136.09  # g/mol
    mass_Na2HPO4_2H2O = 1.00  # g
    mw_Na2HPO4_2H2O = 177.99  # g/mol
    Ka2 = 6.2e-8
    Ka3 = 1.8e-12

    # --- Options from the Question ---
    # Note: These are the options as listed in the final user prompt, not the shuffled lists from some agents.
    options = {
        'A': 3.97e-7,
        'B': 5.48e-7,
        'C': 2.81e-7,
        'D': 6.24e-7
    }
    
    # --- The Final Answer to be Checked ---
    final_answer_letter = 'D'

    # --- Step 1: Calculate moles of each component ---
    # Moles of H2PO4- from KH2PO4
    moles_H2PO4_minus = mass_KH2PO4 / mw_KH2PO4
    
    # Moles of HPO4^2- from Na2HPO4.2H2O
    moles_HPO4_2minus = mass_Na2HPO4_2H2O / mw_Na2HPO4_2H2O

    # --- Step 2: Calculate initial concentrations ---
    # Convert volume from cm^3 to Liters
    volume_L = volume_cm3 / 1000.0
    
    # Concentration of H2PO4-
    conc_H2PO4_minus = moles_H2PO4_minus / volume_L
    
    # Concentration of HPO4^2-
    conc_HPO4_2minus = moles_HPO4_2minus / volume_L

    # --- Step 3: Calculate the hydrogen ion concentration [H+] ---
    # The buffer system is H2PO4-/HPO4^2-, governed by Ka2.
    # Ka2 = [H+][HPO4^2-] / [H2PO4-]
    # Rearranging gives: [H+] = Ka2 * [H2PO4-] / [HPO4^2-]
    conc_H_plus = Ka2 * (conc_H2PO4_minus / conc_HPO4_2minus)

    # --- Step 4: Calculate the orthophosphate ion concentration [PO4^3-] ---
    # The formation of PO4^3- is governed by Ka3.
    # Ka3 = [H+][PO4^3-] / [HPO4^2-]
    # Rearranging gives: [PO4^3-] = Ka3 * [HPO4^2-] / [H+]
    calculated_conc_PO4_3minus = Ka3 * (conc_HPO4_2minus / conc_H_plus)

    # --- Step 5: Verify the correctness of the final answer ---
    # Get the value corresponding to the final answer's letter
    final_answer_value = options.get(final_answer_letter)

    if final_answer_value is None:
        return f"Invalid answer option '{final_answer_letter}'. Valid options are A, B, C, D."

    # Compare the calculated value with the final answer's value using a relative tolerance.
    # A 1% tolerance is reasonable given the options have 3 significant figures.
    if math.isclose(calculated_conc_PO4_3minus, final_answer_value, rel_tol=0.01):
        return "Correct"
    else:
        # If incorrect, find the correct option and provide a detailed reason.
        correct_letter = 'Unknown'
        for letter, value in options.items():
            if math.isclose(calculated_conc_PO4_3minus, value, rel_tol=0.01):
                correct_letter = letter
                break
        
        reason = (
            f"Incorrect. The provided answer is {final_answer_letter} ({final_answer_value:.2e} M), "
            f"but the calculated concentration of orthophosphate ions is {calculated_conc_PO4_3minus:.3e} M. "
            f"This value correctly matches option {correct_letter} ({options.get(correct_letter):.2e} M).\n\n"
            f"Calculation Breakdown:\n"
            f"1. Moles H₂PO₄⁻ = {moles_H2PO4_minus:.4e} mol\n"
            f"2. Moles HPO₄²⁻ = {moles_HPO4_2minus:.4e} mol\n"
            f"3. [H₂PO₄⁻] = {conc_H2PO4_minus:.4e} M\n"
            f"4. [HPO₄²⁻] = {conc_HPO4_2minus:.4e} M\n"
            f"5. [H⁺] = {conc_H_plus:.4e} M (calculated using Ka₂)\n"
            f"6. [PO₄³⁻] = {calculated_conc_PO4_3minus:.4e} M (calculated using Ka₃)"
        )
        return reason

# Execute the checking function and print the result.
result = check_chemistry_answer()
print(result)