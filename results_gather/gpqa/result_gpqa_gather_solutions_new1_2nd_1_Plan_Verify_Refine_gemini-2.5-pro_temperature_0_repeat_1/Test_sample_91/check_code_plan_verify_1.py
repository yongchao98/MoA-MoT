import math

def check_enthalpy_calculation():
    """
    This function checks the calculation for the enthalpy of formation of
    (CH3)2C=CH-CH2-CH(CH3)-CH2-CH=C(CH3)2 based on the provided data.
    It verifies the steps of the calculation and compares the result with the given options.
    """
    # --- Given Data from the Question ---
    h_atom_C = 1000  # Enthalpy of atomization of carbon (kJ/mol)
    be_H_H = 100     # Bond energy of H-H (kJ/mol)
    be_C_C = 200     # Bond energy of C-C (kJ/mol)
    be_C_eq_C = 300  # Bond energy of C=C (kJ/mol)
    be_C_H = 400     # Bond energy of C-H (kJ/mol)

    # --- Options from the Question ---
    options = {
        "A": "67.4 kJ/mol",
        "B": "1900 kJ/g",
        "C": "11.44 kJ/g",
        "D": "11200 kJ/mol"
    }
    
    # The final answer provided by the LLM to be checked
    llm_final_answer = "C"

    # --- Step 1: Determine Molecular Formula and Bond Counts ---
    # Molecule: (CH3)2C=CH-CH2-CH(CH3)-CH2-CH=C(CH3)2
    # Molecular Formula is C12H22
    num_C_atoms = 12
    num_H_atoms = 22
    
    # Bond counts derived from the structure
    num_C_H_bonds = 22
    num_C_C_bonds = 9
    num_C_eq_C_bonds = 2

    # --- Step 2: Calculate the Enthalpy of Atomization of Reactants ---
    # The formation reaction is: 12 C(s) + 11 H2(g) -> C12H22(g)
    # Energy input = (12 * ΔH_atom(C)) + (11 * BE(H-H))
    num_H2_moles = num_H_atoms / 2
    enthalpy_atomization_reactants = (num_C_atoms * h_atom_C) + (num_H2_moles * be_H_H)

    # --- Step 3: Calculate the Enthalpy of Atomization of the Product ---
    # This is the sum of all bond energies in the product molecule.
    enthalpy_atomization_product = (num_C_C_bonds * be_C_C) + \
                                    (num_C_eq_C_bonds * be_C_eq_C) + \
                                    (num_C_H_bonds * be_C_H)

    # --- Step 4: Calculate the Enthalpy of Formation (ΔHf) in kJ/mol ---
    # ΔHf = (Energy to atomize reactants) - (Energy released forming product)
    delta_hf_mol = enthalpy_atomization_reactants - enthalpy_atomization_product

    # --- Step 5: Convert to kJ/g to Compare with Options ---
    # Molar mass of C12H22 (using integer masses as implied by the problem's simple numbers)
    molar_mass_C = 12  # g/mol
    molar_mass_H = 1   # g/mol
    molar_mass_C12H22 = (num_C_atoms * molar_mass_C) + (num_H_atoms * molar_mass_H)
    
    delta_hf_g = delta_hf_mol / molar_mass_C12H22

    # --- Verification ---
    # Check if the LLM's chosen answer is correct
    correct_value_str = options.get(llm_final_answer)
    if not correct_value_str:
        return f"Incorrect. The final answer '{llm_final_answer}' is not a valid option."
        
    correct_value = float(correct_value_str.split()[0])
    correct_unit = correct_value_str.split()[1]

    # Check the value and unit
    if correct_unit == "kJ/g":
        calculated_value = delta_hf_g
    elif correct_unit == "kJ/mol":
        calculated_value = delta_hf_mol
    else:
        return f"Incorrect. The unit '{correct_unit}' in the chosen answer is not recognized."

    # Use a small tolerance for floating point comparison
    if not math.isclose(calculated_value, correct_value, rel_tol=1e-2):
        return f"Incorrect. The final answer is {llm_final_answer} ({correct_value_str}), but the calculated value is {calculated_value:.2f} {correct_unit}. " \
               f"Calculation details: Reactant atomization energy = {enthalpy_atomization_reactants} kJ. " \
               f"Product bond energy = {enthalpy_atomization_product} kJ/mol. " \
               f"Enthalpy of formation = {delta_hf_mol} kJ/mol. " \
               f"Molar mass = {molar_mass_C12H22} g/mol."

    # Check for common distractors to ensure the reasoning is sound
    # Option D is the product atomization energy
    if math.isclose(float(options["D"].split()[0]), enthalpy_atomization_product):
        pass # This is expected
    else:
        return f"Logic check failed: The value of option D ({options['D']}) should be the product atomization energy, but the calculated value is {enthalpy_atomization_product} kJ/mol."

    # Option B has the correct value for kJ/mol but the wrong unit
    if math.isclose(float(options["B"].split()[0]), delta_hf_mol):
        pass # This is expected
    else:
        return f"Logic check failed: The numerical value of option B ({options['B']}) should match the enthalpy of formation in kJ/mol, but the calculated value is {delta_hf_mol} kJ/mol."

    return "Correct"

# Run the check
result = check_enthalpy_calculation()
print(result)