import math

def check_answer():
    """
    This function checks the correctness of the provided LLM answer by recalculating the enthalpy of formation from the given data.
    """
    # --- 1. Define problem constants from the question ---
    dHa_C = 1000      # Enthalpy of atomization of Carbon (kJ/mol)
    BE_HH = 100       # Bond energy of H-H (kJ/mol)
    BE_CC = 200       # Bond energy of C-C (kJ/mol)
    BE_CeqC = 300     # Bond energy of C=C (kJ/mol)
    BE_CH = 400       # Bond energy of C-H (kJ/mol)
    
    llm_provided_answer = "C"

    # --- 2. Determine molecular properties ---
    # Molecule: (CH3)2C=CH-CH2-CH(CH3)-CH2-CH=C(CH3)2
    # Formula is C12H22. Formation reaction: 12 C(s) + 11 H2(g) -> C12H22(g)
    num_C_atoms = 12
    num_H2_molecules = 11
    
    # Bond counts in the product molecule
    num_CC_bonds = 9
    num_CeqC_bonds = 2
    num_CH_bonds = 22

    # --- 3. Perform the enthalpy calculation ---
    # ΔHf° = [Σ (Enthalpy of atomization of reactants)] - [Σ (Bond energies of product)]
    
    # Enthalpy of atomization of reactants: 12 C(s) + 11 H2(g) -> 12 C(g) + 22 H(g)
    dHa_reactants = (num_C_atoms * dHa_C) + (num_H2_molecules * BE_HH)
    
    # Sum of bond energies in the product molecule (its enthalpy of atomization)
    dHa_product = (num_CC_bonds * BE_CC) + (num_CeqC_bonds * BE_CeqC) + (num_CH_bonds * BE_CH)
    
    # Enthalpy of formation in kJ/mol
    dHf_kJ_per_mol = dHa_reactants - dHa_product

    # --- 4. Convert units to check all options ---
    # Molar mass of C12H22 (using integer masses consistent with problem data)
    molar_mass = (12 * 12.0) + (22 * 1.0)
    
    # Enthalpy of formation in kJ/g
    dHf_kJ_per_g = dHf_kJ_per_mol / molar_mass

    # --- 5. Verify the LLM's answer against the calculated result ---
    options = {
        "A": {"value": 67.4, "unit": "kJ/mol"},
        "B": {"value": 11200, "unit": "kJ/mol"},
        "C": {"value": 11.44, "unit": "kJ/g"},
        "D": {"value": 1900, "unit": "kJ/g"}
    }

    correct_option = None
    # A small tolerance is used for floating point comparisons
    tolerance = 1e-2 

    # Check which option matches the calculation
    if math.isclose(dHf_kJ_per_g, options["C"]["value"], rel_tol=tolerance):
        correct_option = "C"
    elif math.isclose(dHf_kJ_per_mol, options["D"]["value"], rel_tol=tolerance):
        # This case checks if the unit was mistaken
        correct_option = "D" 
    elif math.isclose(dHf_kJ_per_mol, options["A"]["value"], rel_tol=tolerance):
        correct_option = "A"
    elif math.isclose(dHf_kJ_per_mol, options["B"]["value"], rel_tol=tolerance):
        correct_option = "B"

    if correct_option == llm_provided_answer:
        return "Correct"
    else:
        # Provide a detailed reason for the discrepancy
        reason = (
            f"Incorrect. The provided answer is {llm_provided_answer}, but the calculation points to option {correct_option}.\n"
            f"Calculation Steps:\n"
            f"1. Reactant atomization energy = (12 * {dHa_C}) + (11 * {BE_HH}) = {dHa_reactants} kJ.\n"
            f"2. Product bond energy sum = (9 * {BE_CC}) + (2 * {BE_CeqC}) + (22 * {BE_CH}) = {dHa_product} kJ.\n"
            f"3. Enthalpy of formation (kJ/mol) = {dHa_reactants} - {dHa_product} = {dHf_kJ_per_mol} kJ/mol.\n"
            f"4. Molar mass of C12H22 = {molar_mass} g/mol.\n"
            f"5. Enthalpy of formation (kJ/g) = {dHf_kJ_per_mol} / {molar_mass} = {dHf_kJ_per_g:.2f} kJ/g.\n"
            f"This calculated value of {dHf_kJ_per_g:.2f} kJ/g matches option C ({options['C']['value']} kJ/g)."
        )
        return reason

# Execute the check and print the result
result = check_answer()
print(result)