import math

def check_correctness():
    """
    This function checks the correctness of the calculated enthalpy of formation for the given molecule.
    It follows the standard procedure using bond energies and enthalpy of atomization.
    """
    
    # --- Given Data from the Question ---
    enthalpy_atomization_C = 1000  # kJ/mol
    bond_energy_H_H = 100         # kJ/mol
    bond_energy_C_C = 200         # kJ/mol
    bond_energy_C_eq_C = 300      # kJ/mol
    bond_energy_C_H = 400         # kJ/mol
    
    # --- Analysis of the Molecule: (CH3)2C=CH-CH2-CH(CH3)-CH2-CH=C(CH3)2 ---
    
    # Step 1: Determine the molecular formula (C_n H_m)
    # C atoms: 2+1+1+1+(1+1)+1+1+(1+2) = 12
    # H atoms: 6+1+2+(1+3)+2+1+6 = 22
    num_C = 12
    num_H = 22
    
    # Step 2: Count the number of each type of bond
    # C-H bonds = number of H atoms
    num_C_H_bonds = 22
    # C=C bonds are explicitly shown
    num_C_eq_C_bonds = 2
    # For an acyclic molecule, total C-C linkages = num_C - 1.
    # C-C single bonds = total linkages - double bonds
    num_C_C_bonds = num_C - 1 - num_C_eq_C_bonds
    
    # --- Calculation ---
    
    # Step 3: Calculate the enthalpy of atomization of reactants
    # The formation reaction is: 12 C(s) + 11 H2(g) -> C12H22(g)
    # Energy to atomize reactants = (12 * ΔH_atom(C)) + (11 * BE(H-H))
    num_H2_moles = num_H / 2
    enthalpy_atomization_reactants = (num_C * enthalpy_atomization_C) + (num_H2_moles * bond_energy_H_H)
    
    # Step 4: Calculate the total bond energy of the product
    # This is the energy released when forming the molecule from gaseous atoms.
    total_bond_energy_product = (num_C_C_bonds * bond_energy_C_C) + \
                                (num_C_eq_C_bonds * bond_energy_C_eq_C) + \
                                (num_C_H_bonds * bond_energy_C_H)
                                
    # Step 5: Calculate the enthalpy of formation (ΔH_f) in kJ/mol
    # ΔH_f = (Energy to break reactant bonds) - (Energy released forming product bonds)
    delta_H_f_mol = enthalpy_atomization_reactants - total_bond_energy_product
    
    # Step 6: Convert to kJ/g
    # Molar mass is calculated using integer atomic weights as implied by the problem context.
    molar_mass_C = 12  # g/mol
    molar_mass_H = 1   # g/mol
    molar_mass_compound = (num_C * molar_mass_C) + (num_H * molar_mass_H)
    
    if molar_mass_compound == 0:
        return "Error: Molar mass is zero, cannot divide."
        
    delta_H_f_gram = delta_H_f_mol / molar_mass_compound
    
    # --- Verification ---
    
    # The provided answer is 'C', which corresponds to option C) 11.44 kJ/g
    expected_answer_value = 11.44
    
    # Check if the calculated value matches the expected answer within a small tolerance.
    if math.isclose(delta_H_f_gram, expected_answer_value, rel_tol=1e-3):
        return "Correct"
    else:
        # Provide a detailed reason for the incorrectness.
        # Check for common distractors to provide a more specific error message.
        if math.isclose(total_bond_energy_product, 11200):
             return (f"Incorrect. The final answer is wrong because it seems to be the enthalpy of atomization of the product, not the enthalpy of formation. "
                     f"The calculated enthalpy of atomization of the product is {total_bond_energy_product} kJ/mol, which matches option A.")
        
        if math.isclose(delta_H_f_mol, 1900):
             return (f"Incorrect. The final answer is wrong because of a unit mismatch. "
                     f"The calculated enthalpy of formation is {delta_H_f_mol} kJ/mol. "
                     f"Option B has the value 1900 but with units of kJ/g. The correct value in kJ/g is {delta_H_f_gram:.2f} kJ/g.")

        return (f"Incorrect. The calculated enthalpy of formation is {delta_H_f_gram:.2f} kJ/g, which does not match the expected answer of {expected_answer_value} kJ/g.\n"
                f"Calculation details:\n"
                f"- Reactant atomization energy: {enthalpy_atomization_reactants} kJ\n"
                f"- Product bond energy: {total_bond_energy_product} kJ/mol\n"
                f"- Enthalpy of formation (mol): {delta_H_f_mol} kJ/mol\n"
                f"- Molar mass: {molar_mass_compound} g/mol\n"
                f"- Enthalpy of formation (gram): {delta_H_f_gram:.4f} kJ/g")

# Execute the check and print the result
result = check_correctness()
print(result)