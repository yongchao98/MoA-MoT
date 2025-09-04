import math

def check_enthalpy_of_formation():
    """
    This function checks the correctness of the given answer for the enthalpy of formation calculation.
    It recalculates the value based on the provided data and compares it to the answer.
    """
    # --- Given Constants from the Question ---
    enthalpy_atom_C = 1000  # kJ/mol, for C(s) -> C(g)
    be_H_H = 100          # kJ/mol, for H2(g) -> 2H(g)
    be_C_C = 200          # kJ/mol
    be_C_eq_C = 300       # kJ/mol (C=C double bond)
    be_C_H = 400          # kJ/mol

    # --- Step 1: Analyze the Molecule ---
    # Molecule: (CH3)2C=CH-CH2-CH(CH3)-CH2-CH=C(CH3)2
    # Molecular Formula: C12H22
    num_C_atoms = 12
    num_H_atoms = 22

    # Bond Counts:
    # - There are two C=C double bonds.
    # - The total number of bonds in the carbon skeleton of an acyclic molecule is (num_C_atoms - 1).
    # - Total C-C and C=C bonds = 12 - 1 = 11.
    # - Since there are 2 C=C bonds, the number of C-C single bonds is 11 - 2 = 9.
    # - The number of C-H bonds is equal to the number of hydrogen atoms.
    num_C_H_bonds = 22
    num_C_eq_C_bonds = 2
    num_C_C_bonds = 9

    # --- Step 2: Calculate Enthalpy of Atomization of Reactants ---
    # The formation reaction from elements in their standard state is:
    # 12 C(s) + 11 H2(g) -> C12H22(g)
    # Energy to atomize reactants = (Energy to atomize 12 C) + (Energy to atomize 11 H2)
    num_H2_moles = num_H_atoms / 2
    enthalpy_atomization_reactants = (num_C_atoms * enthalpy_atom_C) + (num_H2_moles * be_H_H)

    # --- Step 3: Calculate Total Bond Energy of the Product ---
    # This is the energy released when gaseous atoms form the molecule.
    # It is also known as the enthalpy of atomization of the product.
    total_bond_energy_product = (num_C_C_bonds * be_C_C) + \
                                (num_C_eq_C_bonds * be_C_eq_C) + \
                                (num_C_H_bonds * be_C_H)

    # --- Step 4: Calculate Enthalpy of Formation (ΔH_f) in kJ/mol ---
    # ΔH_f = (Enthalpy of atomization of reactants) - (Enthalpy of atomization of product)
    delta_H_f_mol = enthalpy_atomization_reactants - total_bond_energy_product

    # --- Step 5: Convert to kJ/g for Comparison ---
    # The given answer 'B' is in kJ/g.
    # We will use integer atomic masses (C=12, H=1) as is common in such problems to match the answer's precision.
    molar_mass_C = 12
    molar_mass_H = 1
    total_molar_mass = (num_C_atoms * molar_mass_C) + (num_H_atoms * molar_mass_H)
    
    if total_molar_mass == 0:
        return "Error: Molar mass cannot be zero."

    calculated_delta_H_f_gram = delta_H_f_mol / total_molar_mass

    # --- Step 6: Verify the Answer ---
    # The provided answer is B, which corresponds to 11.44 kJ/g.
    answer_b_value = 11.44
    
    # Check if the calculated value is close to the answer's value (allowing for small rounding differences).
    if math.isclose(calculated_delta_H_f_gram, answer_b_value, rel_tol=1e-2):
        return "Correct"
    else:
        # If the answer is incorrect, provide a detailed reason.
        reason = (
            f"The provided answer is incorrect. Here is the step-by-step calculation:\n"
            f"1.  **Reactant Atomization Enthalpy**: For 12 C(s) + 11 H2(g) -> 12 C(g) + 22 H(g)\n"
            f"    Energy = (12 * {enthalpy_atom_C}) + (11 * {be_H_H}) = {enthalpy_atomization_reactants} kJ\n"
            f"2.  **Product Bond Energy**: For C12H22, with 9 C-C, 2 C=C, and 22 C-H bonds.\n"
            f"    Energy = (9 * {be_C_C}) + (2 * {be_C_eq_C}) + (22 * {be_C_H}) = {total_bond_energy_product} kJ\n"
            f"3.  **Enthalpy of Formation (kJ/mol)**: ΔH_f = {enthalpy_atomization_reactants} - {total_bond_energy_product} = {delta_H_f_mol} kJ/mol\n"
            f"4.  **Molar Mass**: (12 * 12) + (22 * 1) = {total_molar_mass} g/mol\n"
            f"5.  **Enthalpy of Formation (kJ/g)**: {delta_H_f_mol} / {total_molar_mass} = {calculated_delta_H_f_gram:.4f} kJ/g\n"
            f"The calculated value is {calculated_delta_H_f_gram:.4f} kJ/g, which does not match the value from answer B (11.44 kJ/g)."
        )
        return reason

# Execute the check and print the result
result = check_enthalpy_of_formation()
print(result)