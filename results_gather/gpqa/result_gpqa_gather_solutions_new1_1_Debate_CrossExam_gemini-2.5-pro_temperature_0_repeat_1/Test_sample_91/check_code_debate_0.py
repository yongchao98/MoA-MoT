import math

def check_answer():
    """
    This function checks the correctness of the calculated enthalpy of formation.
    """
    # Given values from the question
    H_atom_C = 1000  # Enthalpy of atomization of carbon, kJ/mol
    BE_HH = 100      # Bond energy of H-H, kJ/mol
    BE_CC = 200      # Bond energy of C-C, kJ/mol
    BE_C_dbl_C = 300 # Bond energy of C=C, kJ/mol
    BE_CH = 400      # Bond energy of C-H, kJ/mol

    # Step 1: Determine Molecular Formula and Bond Count for (CH3)2C=CH-CH2-CH(CH3)-CH2-CH=C(CH3)2
    # By analyzing the structure:
    num_C = 12
    num_H = 22
    num_CH_bonds = 22
    num_CC_single_bonds = 9
    num_CC_double_bonds = 2

    # Step 2: Calculate the enthalpy of atomization for the reactants
    # The formation reaction is: 12 C(s) + 11 H2(g) -> C12H22(g)
    # Energy required to form gaseous atoms from elements in their standard state
    moles_of_H2 = num_H / 2
    H_atom_reactants = (num_C * H_atom_C) + (moles_of_H2 * BE_HH)

    # Step 3: Calculate the sum of bond energies for the product molecule (C12H22)
    # This is the energy released when forming the molecule from its gaseous atoms.
    # It is also known as the enthalpy of atomization of the product.
    H_atom_product = (num_CH_bonds * BE_CH) + (num_CC_single_bonds * BE_CC) + (num_CC_double_bonds * BE_C_dbl_C)

    # Step 4: Calculate the enthalpy of formation (in kJ/mol)
    # ΔHf = (Energy to atomize reactants) - (Energy released forming product)
    delta_H_f_mol = H_atom_reactants - H_atom_product

    # Step 5: Convert the enthalpy of formation to kJ/g
    # Calculate the molar mass of C12H22 using integer masses as is common for such problems.
    molar_mass_C = 12  # g/mol
    molar_mass_H = 1   # g/mol
    molar_mass_compound = (num_C * molar_mass_C) + (num_H * molar_mass_H)
    
    if molar_mass_compound == 0:
        return "Error: Molar mass is zero, cannot divide."
        
    delta_H_f_gram = delta_H_f_mol / molar_mass_compound

    # The final answer provided is 'C', which corresponds to 11.44 kJ/g.
    # Let's check if our calculation matches this value.
    expected_answer_C = 11.44
    
    # Check for common distractors
    option_B_value = 1900 # This is the value for ΔHf in kJ/g
    option_D_value = 11200 # This is the value for the product's enthalpy of atomization in kJ/mol

    if not math.isclose(delta_H_f_gram, expected_answer_C, rel_tol=1e-3):
        reason = f"The calculated enthalpy of formation is {delta_H_f_gram:.2f} kJ/g, which does not match the value in option C (11.44 kJ/g).\n"
        reason += "Let's re-check the steps:\n"
        reason += f"1. Reactant atomization energy: (12 * {H_atom_C}) + (11 * {BE_HH}) = {H_atom_reactants} kJ\n"
        reason += f"2. Product bond energy sum: (9 * {BE_CC}) + (2 * {BE_C_dbl_C}) + (22 * {BE_CH}) = {H_atom_product} kJ\n"
        reason += f"3. Enthalpy of formation (kJ/mol): {H_atom_reactants} - {H_atom_product} = {delta_H_f_mol} kJ/mol\n"
        reason += f"4. Molar Mass: {molar_mass_compound} g/mol\n"
        reason += f"5. Enthalpy of formation (kJ/g): {delta_H_f_mol} / {molar_mass_compound} = {delta_H_f_gram:.4f} kJ/g\n"
        reason += "The calculation is correct. The discrepancy might be due to a typo in the provided answer or options."
        return reason

    # The calculation matches the provided answer.
    return "Correct"

# Run the check
result = check_answer()
print(result)