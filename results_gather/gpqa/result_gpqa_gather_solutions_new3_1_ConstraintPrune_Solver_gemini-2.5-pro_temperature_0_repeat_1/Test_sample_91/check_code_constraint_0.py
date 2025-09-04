import math

def check_enthalpy_of_formation():
    """
    This function calculates the enthalpy of formation for the given molecule and checks
    if it matches the provided answer 'C'.

    The calculation follows these steps:
    1.  Define the molecular structure and count all atoms and bonds.
    2.  Use the provided thermodynamic data.
    3.  Calculate the total enthalpy of atomization for the reactants (elements in standard state).
    4.  Calculate the total bond energy (enthalpy of atomization) for the product molecule.
    5.  Calculate the enthalpy of formation (ΔH_f) in kJ/mol.
    6.  Convert ΔH_f to kJ/g using the molecule's molar mass.
    7.  Compare the calculated value to the value in option C.
    """
    # --- Given Thermodynamic Data ---
    enthalpy_atomization_C = 1000  # kJ/mol
    bond_energy_HH = 100           # kJ/mol
    bond_energy_CC = 200           # kJ/mol
    bond_energy_C_eq_C = 300       # kJ/mol
    bond_energy_CH = 400           # kJ/mol

    # --- Step 1: Molecular Formula and Bond Counts ---
    # Molecule: (CH3)2C=CH-CH2-CH(CH3)-CH2-CH=C(CH3)2
    num_C = 12
    num_H = 22
    num_CH_bonds = 22
    num_C_eq_C_bonds = 2
    num_CC_bonds = 9  # Verified by counting: 2+1+1+1+1+1+2 = 9

    # --- Step 2: Calculate Enthalpy of Atomization of Reactants ---
    # Formation reaction: 12 C(s) + 11 H2(g) -> C12H22(g)
    # Energy to form gaseous atoms: 12 C(g) + 22 H(g)
    enthalpy_atom_reactants = (num_C * enthalpy_atomization_C) + ((num_H / 2) * bond_energy_HH)
    
    # --- Step 3: Calculate Enthalpy of Atomization of Product ---
    # This is the sum of all bond energies in the molecule.
    enthalpy_atom_product = (num_CC_bonds * bond_energy_CC) + \
                            (num_C_eq_C_bonds * bond_energy_C_eq_C) + \
                            (num_CH_bonds * bond_energy_CH)

    # --- Step 4: Calculate Enthalpy of Formation (kJ/mol) ---
    delta_H_f_mol = enthalpy_atom_reactants - enthalpy_atom_product

    # --- Step 5: Convert to kJ/g ---
    # Molar Mass of C12H22 (using integer masses C=12, H=1 as is standard for such problems)
    molar_mass = (num_C * 12) + (num_H * 1)
    
    if molar_mass == 0:
        return "Error: Molar mass cannot be zero."
        
    delta_H_f_gram = delta_H_f_mol / molar_mass

    # --- Step 6: Check against the provided answer ---
    # The provided answer is 'C', which corresponds to 11.44 kJ/g.
    expected_value_from_option_C = 11.44
    
    # Check if the calculated value is close to the expected value from option C.
    # A relative tolerance of 0.5% is used to account for potential rounding in the problem's options.
    if math.isclose(delta_H_f_gram, expected_value_from_option_C, rel_tol=0.005):
        return "Correct"
    else:
        reason = (
            f"The provided answer 'C' (11.44 kJ/g) is incorrect.\n"
            f"The calculation is as follows:\n"
            f"1. The molecular formula is C12H22, with 9 C-C, 2 C=C, and 22 C-H bonds.\n"
            f"2. Enthalpy of atomization of reactants (12 C(s) + 11 H2(g)) = {enthalpy_atom_reactants} kJ.\n"
            f"3. Enthalpy of atomization of the product (sum of bond energies) = {enthalpy_atom_product} kJ/mol.\n"
            f"4. Enthalpy of formation (ΔH_f) = {enthalpy_atom_reactants} - {enthalpy_atom_product} = {delta_H_f_mol} kJ/mol.\n"
            f"5. Molar mass = {molar_mass} g/mol.\n"
            f"6. Final calculated ΔH_f = {delta_H_f_mol} / {molar_mass} = {delta_H_f_gram:.4f} kJ/g.\n"
            f"This calculated value does not match the value in option C."
        )
        return reason

# Execute the check and print the result
result = check_enthalpy_of_formation()
print(result)