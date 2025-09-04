import re
import math

def check_enthalpy_calculation():
    """
    This function checks the correctness of the calculated enthalpy of formation for the given molecule.
    It performs the calculation from scratch based on the provided data and compares it to the selected answer.
    """
    # --- Given Data from the Question ---
    H_atom_C = 1000      # Enthalpy of atomization of C(s) -> C(g) in kJ/mol
    BE_HH = 100          # Bond energy of H-H in kJ/mol
    BE_CC = 200          # Bond energy of C-C in kJ/mol
    BE_C_dbl_C = 300     # Bond energy of C=C in kJ/mol
    BE_CH = 400          # Bond energy of C-H in kJ/mol

    # --- Step 1: Analyze the Molecule ---
    # Molecule: (CH3)2C=CH-CH2-CH(CH3)-CH2-CH=C(CH3)2
    # Molecular Formula: C12H22
    num_C = 12
    num_H = 22
    # Bond Counts:
    # C-H bonds = number of H atoms = 22
    # C=C bonds = 2 (explicitly shown)
    # C-C bonds = (num_C - 1) - num_C=C = (12 - 1) - 2 = 9
    num_CH_bonds = 22
    num_CC_bonds = 9
    num_C_dbl_C_bonds = 2

    # --- Step 2: Calculate Enthalpy of Atomization of Reactants ---
    # The formation reaction is: 12 C(s) + 11 H2(g) -> C12H22(g)
    # Energy to atomize reactants = 12 * ΔHa(C) + 11 * BE(H-H)
    H_atom_reactants = (num_C * H_atom_C) + ((num_H / 2) * BE_HH)

    # --- Step 3: Calculate Enthalpy of Atomization of Product ---
    # This is the sum of all bond energies in the product molecule.
    H_atom_product = (num_CH_bonds * BE_CH) + (num_CC_bonds * BE_CC) + (num_C_dbl_C_bonds * BE_C_dbl_C)

    # --- Step 4: Calculate Enthalpy of Formation (ΔHf) in kJ/mol ---
    # ΔHf = (Energy to atomize reactants) - (Energy released forming product)
    delta_H_f_mol = H_atom_reactants - H_atom_product

    # --- Step 5: Convert to kJ/g if necessary ---
    # Molar Mass of C12H22 (using integer masses for C=12, H=1)
    molar_mass_compound = (num_C * 12) + (num_H * 1)
    # Enthalpy of formation in kJ/g
    delta_H_f_gram = delta_H_f_mol / molar_mass_compound

    # --- Step 6: Check the LLM's Answer ---
    llm_answer_text = "<<<D>>>"
    
    # Define the options from the question
    options = {
        "A": {"value": 67.4, "unit": "kJ/mol"},
        "B": {"value": 1900, "unit": "kJ/g"},
        "C": {"value": 11200, "unit": "kJ/mol"},
        "D": {"value": 11.44, "unit": "kJ/g"}
    }

    # Extract the letter from the LLM's answer
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Invalid answer format. Expected <<<A>>>, <<<B>>>, <<<C>>>, or <<<D>>>."
    
    selected_option_letter = match.group(1)
    selected_option = options[selected_option_letter]
    
    # --- Final Verification ---
    # Check for common pitfalls first
    if selected_option_letter == 'C':
        if math.isclose(H_atom_product, selected_option["value"]):
            return f"Incorrect. The selected answer C ({selected_option['value']} kJ/mol) represents the total bond energy (enthalpy of atomization) of the product, not the enthalpy of formation as requested by the question."
    
    if selected_option_letter == 'B':
        if math.isclose(delta_H_f_mol, selected_option["value"]):
            return f"Incorrect. The selected answer B has the correct numerical value for the enthalpy of formation in kJ/mol ({delta_H_f_mol} kJ/mol), but the units are wrong (kJ/g instead of kJ/mol)."

    # Check if the selected option is the correct final answer
    if selected_option_letter == 'D':
        if math.isclose(delta_H_f_gram, selected_option["value"], rel_tol=1e-2):
            return "Correct"
        else:
            # This case handles if the letter is D but the value is somehow wrong
            return f"Incorrect. The selected option is D, but the calculated value is {delta_H_f_gram:.2f} kJ/g, which does not match the option's value of {selected_option['value']} kJ/g."

    # If any other option is chosen, it's numerically incorrect.
    return (f"Incorrect. The selected answer {selected_option_letter} ({selected_option['value']} {selected_option['unit']}) is not correct.\n"
            f"The calculated enthalpy of formation is {delta_H_f_mol} kJ/mol, which is equivalent to {delta_H_f_gram:.2f} kJ/g.\n"
            f"The correct choice is D (11.44 kJ/g).")

# Run the check
result = check_enthalpy_calculation()
print(result)