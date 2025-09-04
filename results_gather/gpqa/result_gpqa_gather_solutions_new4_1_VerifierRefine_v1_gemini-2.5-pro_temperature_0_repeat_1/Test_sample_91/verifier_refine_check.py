import math

def check_answer():
    """
    This function checks the correctness of the calculated enthalpy of formation for the given molecule.
    """
    # --- Given Data from the Question ---
    enthalpy_atomization_C = 1000  # kJ/mol
    bond_energy_HH = 100           # kJ/mol
    bond_energy_CC = 200           # kJ/mol
    bond_energy_C_dbl_C = 300      # kJ/mol
    bond_energy_CH = 400           # kJ/mol

    # --- Step 1: Determine Molecular Formula and Bond Counts ---
    # Molecule: (CH3)2C=CH-CH2-CH(CH3)-CH2-CH=C(CH3)2
    # By inspection, the molecular formula is C12H22.
    num_C = 12
    num_H = 22
    
    # Counting the bonds:
    # C-H bonds = number of H atoms
    num_CH_bonds = 22
    # C=C bonds are explicitly shown
    num_C_dbl_C_bonds = 2
    # C-C single bonds: For an acyclic molecule with 12 carbons, there are 11 total C-C links.
    # Since 2 are double bonds, the remaining 9 must be single bonds.
    num_CC_bonds = 9

    # --- Step 2: Calculate Enthalpy of Atomization of Reactants ---
    # The formation reaction is: 12 C(s) + 11 H2(g) -> C12H22(g)
    # Energy to convert reactants to gaseous atoms: 12 C(g) + 22 H(g)
    # We need 11 moles of H2 to get 22 moles of H atoms.
    enthalpy_atomization_reactants = (num_C * enthalpy_atomization_C) + ((num_H / 2) * bond_energy_HH)
    
    # --- Step 3: Calculate Total Bond Energy of the Product ---
    # This is the energy released when forming one mole of the product from gaseous atoms.
    # It is also the enthalpy of atomization of the product.
    total_bond_energy_product = (num_CC_bonds * bond_energy_CC) + \
                                (num_C_dbl_C_bonds * bond_energy_C_dbl_C) + \
                                (num_CH_bonds * bond_energy_CH)

    # --- Step 4: Calculate Enthalpy of Formation (ΔH_f) in kJ/mol ---
    # ΔH_f = (Enthalpy of Atomization of Reactants) - (Total Bond Energy of Product)
    delta_H_f_mol = enthalpy_atomization_reactants - total_bond_energy_product

    # --- Step 5: Convert to kJ/g to Match the Options ---
    # Calculate molar mass of C12H22 (using integer masses C=12, H=1)
    molar_mass = (num_C * 12) + (num_H * 1)
    
    # Convert ΔH_f from kJ/mol to kJ/g
    if molar_mass == 0:
        return "Error: Molar mass is zero, cannot divide."
    delta_H_f_gram = delta_H_f_mol / molar_mass

    # --- Step 6: Check against the provided answer ---
    # The provided answer is 'D', which corresponds to 11.44 kJ/g.
    expected_answer_value = 11.44
    tolerance = 0.01

    # Check if the calculated value matches the expected answer.
    if not math.isclose(delta_H_f_gram, expected_answer_value, rel_tol=tolerance):
        return (f"Incorrect. The calculated enthalpy of formation is {delta_H_f_mol:.2f} kJ/mol. "
                f"This converts to {delta_H_f_gram:.2f} kJ/g using a molar mass of {molar_mass} g/mol. "
                f"This value does not match the selected answer D (11.44 kJ/g).")

    # Let's also check for common mistakes that match other options.
    # Option C is 11200 kJ/mol, which is the total bond energy of the product.
    if math.isclose(total_bond_energy_product, 11200, rel_tol=tolerance):
        pass # This is a known distractor, but the final answer D is correct.
    
    # Option B is 1900 kJ/g. The value 1900 corresponds to kJ/mol, not kJ/g.
    if math.isclose(delta_H_f_mol, 1900, rel_tol=tolerance):
        if math.isclose(delta_H_f_gram, 1900, rel_tol=tolerance):
            return "Incorrect. The value 1900 corresponds to kJ/mol, but option B incorrectly labels it as kJ/g."

    # If all checks pass and the calculated value matches option D, the answer is correct.
    return "Correct"

# Run the check
result = check_answer()
print(result)