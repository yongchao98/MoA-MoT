import math

def check_enthalpy_calculation():
    """
    This function checks the correctness of the calculated enthalpy of formation for the given molecule.
    It recalculates the value based on the problem's constraints and compares it to the provided answer.
    """
    # Given data from the question
    H_atom_C = 1000      # Enthalpy of atomization of C(s) -> C(g) in kJ/mol
    BE_HH = 100          # Bond energy of H-H in kJ/mol
    BE_CC = 200          # Bond energy of C-C in kJ/mol
    BE_C_dbl_C = 300     # Bond energy of C=C in kJ/mol
    BE_CH = 400          # Bond energy of C-H in kJ/mol

    # Step 1: Determine the molecular formula and count the bonds for (CH3)2C=CH-CH2-CH(CH3)-CH2-CH=C(CH3)2
    # Atom count:
    # (CH3)2C= -> 3C, 6H
    # =CH- -> 1C, 1H
    # -CH2- -> 1C, 2H
    # -CH(CH3)- -> 2C, 4H
    # -CH2- -> 1C, 2H
    # -CH= -> 1C, 1H
    # =C(CH3)2 -> 3C, 6H
    # Total C = 3+1+1+2+1+1+3 = 12
    # Total H = 6+1+2+4+2+1+6 = 22
    num_C = 12
    num_H = 22
    
    # Bond count:
    # C-H bonds = number of H atoms
    num_CH_bonds = 22
    # C=C bonds are explicitly shown
    num_C_dbl_C_bonds = 2
    # C-C single bonds:
    # (CH3)2C= -> 2
    # =C-C- -> 1
    # -C-C( -> 1
    # C-CH3 -> 1
    # )-C-C- -> 1
    # -C-C= -> 1
    # =C(CH3)2 -> 2
    # Total = 2+1+1+1+1+1+2 = 9
    num_CC_bonds = 9

    # Step 2: Calculate the enthalpy of atomization for the reactants: 12 C(s) + 11 H2(g)
    # This is the energy required to form 12 C(g) and 22 H(g).
    # We need 11 moles of H2 to get 22 moles of H.
    num_H2_moles = num_H / 2
    H_atom_reactants = (num_C * H_atom_C) + (num_H2_moles * BE_HH)

    # Step 3: Calculate the sum of bond energies for the product molecule (C12H22)
    # This is the energy released when forming the molecule from its gaseous atoms.
    sum_bond_energies_product = (num_CH_bonds * BE_CH) + (num_CC_bonds * BE_CC) + (num_C_dbl_C_bonds * BE_C_dbl_C)

    # Step 4: Calculate the enthalpy of formation (in kJ/mol)
    # ΔHf° = (Energy to atomize reactants) - (Energy released forming product)
    delta_H_f_mol = H_atom_reactants - sum_bond_energies_product

    # Step 5: Calculate the molar mass of C12H22 (using integer masses as is common in such problems)
    molar_mass_C = 12  # g/mol
    molar_mass_H = 1   # g/mol
    molar_mass_compound = (num_C * molar_mass_C) + (num_H * molar_mass_H)

    # Step 6: Calculate the enthalpy of formation in kJ/g
    if molar_mass_compound == 0:
        return "Error: Molar mass is zero, cannot divide."
    delta_H_f_gram = delta_H_f_mol / molar_mass_compound

    # Step 7: Check the correctness of the final answer 'A'
    # The final answer provided is 'A', which corresponds to 11.44 kJ/g.
    expected_answer_A = 11.44
    
    # Check bond counts
    if num_C != 12 or num_H != 22:
        return f"Incorrect atom count. Calculated C={num_C}, H={num_H}. Expected C=12, H=22."
    if num_CC_bonds != 9 or num_C_dbl_C_bonds != 2 or num_CH_bonds != 22:
        return f"Incorrect bond count. Calculated C-C={num_CC_bonds}, C=C={num_C_dbl_C_bonds}, C-H={num_CH_bonds}. Expected C-C=9, C=C=2, C-H=22."

    # Check intermediate calculations
    if H_atom_reactants != 13100:
        return f"Incorrect reactant atomization energy. Calculated {H_atom_reactants} kJ, expected 13100 kJ."
    if sum_bond_energies_product != 11200:
        return f"Incorrect product bond energy sum. Calculated {sum_bond_energies_product} kJ/mol, expected 11200 kJ/mol."
    if delta_H_f_mol != 1900:
        return f"Incorrect enthalpy of formation in kJ/mol. Calculated {delta_H_f_mol} kJ/mol, expected 1900 kJ/mol."
    
    # Check final value against option A
    if not math.isclose(delta_H_f_gram, expected_answer_A, rel_tol=1e-3):
        return f"Incorrect final answer. The calculated enthalpy of formation is {delta_H_f_gram:.2f} kJ/g, which does not match option A (11.44 kJ/g)."

    return "Correct"

# Run the check
result = check_enthalpy_calculation()
print(result)