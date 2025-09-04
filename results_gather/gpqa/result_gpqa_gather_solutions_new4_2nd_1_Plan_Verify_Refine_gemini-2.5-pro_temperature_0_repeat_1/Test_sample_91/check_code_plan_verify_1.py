import math

def check_enthalpy_calculation():
    """
    This function checks the calculation of the enthalpy of formation for the given molecule.
    """
    # Given data from the question
    h_atom_c = 1000  # kJ/mol
    be_h_h = 100     # kJ/mol
    be_c_c = 200     # kJ/mol
    be_c_eq_c = 300  # kJ/mol
    be_c_h = 400     # kJ/mol

    # Molecule: (CH3)2C=CH-CH2-CH(CH3)-CH2-CH=C(CH3)2
    # 1. Determine molecular formula and bond counts
    num_c = 12
    num_h = 22
    num_c_h_bonds = 22
    num_c_eq_c_bonds = 2
    # For an acyclic molecule, total C-C linkages = num_c - 1
    # num_c_c_bonds = (num_c - 1) - num_c_eq_c_bonds
    num_c_c_bonds = 11 - 2

    # 2. Calculate enthalpy of atomization of reactants
    # Reaction: 12 C(s) + 11 H2(g) -> C12H22(g)
    # We need 12 C(g) and 22 H(g) atoms.
    # Energy for 12 C(s) -> 12 C(g)
    energy_atomize_c = num_c * h_atom_c
    # Energy for 11 H2(g) -> 22 H(g)
    num_h2_moles = num_h / 2
    energy_atomize_h = num_h2_moles * be_h_h
    total_reactant_atomization_energy = energy_atomize_c + energy_atomize_h

    # 3. Calculate total bond energy of the product
    energy_form_c_c = num_c_c_bonds * be_c_c
    energy_form_c_eq_c = num_c_eq_c_bonds * be_c_eq_c
    energy_form_c_h = num_c_h_bonds * be_c_h
    total_product_bond_energy = energy_form_c_c + energy_form_c_eq_c + energy_form_c_h

    # 4. Calculate enthalpy of formation in kJ/mol
    # ΔH_f = Σ(ΔH_atom_reactants) - Σ(BE_products)
    delta_h_f_mol = total_reactant_atomization_energy - total_product_bond_energy

    # 5. Convert to kJ/g
    molar_mass_c = 12  # g/mol
    molar_mass_h = 1   # g/mol
    molar_mass_compound = (num_c * molar_mass_c) + (num_h * molar_mass_h)
    delta_h_f_gram = delta_h_f_mol / molar_mass_compound

    # Options from the original question prompt
    options = {
        "A": 11.44,  # kJ/g
        "B": 67.4,   # kJ/mol
        "C": 1900,   # kJ/g
        "D": 11200   # kJ/mol
    }
    
    # The final answer to check is <<<A>>>
    final_answer_label = "A"
    
    # Check if the calculated value matches the value of the chosen option
    expected_value = options[final_answer_label]
    
    # Check for distractors
    if final_answer_label == "D" and math.isclose(total_product_bond_energy, expected_value, rel_tol=1e-5):
        return f"Incorrect. The answer chose option D, which is {total_product_bond_energy} kJ/mol. This value represents the total bond energy of the product (enthalpy of atomization of the compound), not the enthalpy of formation."
    
    if final_answer_label == "C" and math.isclose(delta_h_f_mol, expected_value, rel_tol=1e-5):
        return f"Incorrect. The answer chose option C, which is {delta_h_f_mol} kJ/g. This uses the correct numerical value for the enthalpy of formation in kJ/mol ({delta_h_f_mol} kJ/mol), but applies the wrong units."

    # Check the final correct answer
    if math.isclose(delta_h_f_gram, expected_value, rel_tol=1e-2):
        return "Correct"
    else:
        return (f"Incorrect. The calculated enthalpy of formation is {delta_h_f_gram:.4f} kJ/g. "
                f"The answer chose option {final_answer_label} with a value of {expected_value} kJ/g. "
                f"The calculation steps are as follows:\n"
                f"1. Reactant atomization energy = (12 * 1000) + (11 * 100) = 13100 kJ\n"
                f"2. Product bond energy = (9 * 200) + (2 * 300) + (22 * 400) = 11200 kJ/mol\n"
                f"3. Enthalpy of formation = 13100 - 11200 = 1900 kJ/mol\n"
                f"4. Molar mass = (12 * 12) + (22 * 1) = 166 g/mol\n"
                f"5. Enthalpy of formation in kJ/g = 1900 / 166 = {delta_h_f_gram:.4f} kJ/g.")

# Run the check
result = check_enthalpy_calculation()
print(result)