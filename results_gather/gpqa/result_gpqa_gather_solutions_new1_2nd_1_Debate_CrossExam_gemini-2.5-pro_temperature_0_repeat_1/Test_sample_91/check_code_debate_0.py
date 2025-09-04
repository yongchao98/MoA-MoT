import math

def check_enthalpy_of_formation():
    """
    This function calculates the enthalpy of formation for the given molecule and
    checks if the provided answer is correct.
    """
    # --- 1. Given Data and Molecular Properties ---
    # Thermodynamic data from the question
    H_atom_C = 1000      # Enthalpy of atomization of Carbon (kJ/mol)
    BE_HH = 100          # Bond energy of H-H (kJ/mol)
    BE_CC_single = 200   # Bond energy of C-C (kJ/mol)
    BE_CC_double = 300   # Bond energy of C=C (kJ/mol)
    BE_CH = 400          # Bond energy of C-H (kJ/mol)

    # Analysis of the molecule: (CH3)2C=CH-CH2-CH(CH3)-CH2-CH=C(CH3)2
    # Molecular Formula: C12H22
    num_C = 12
    num_H = 22
    
    # Bond counts within the molecule
    num_CH_bonds = 22
    num_CC_single_bonds = 9
    num_CC_double_bonds = 2

    # --- 2. Calculate Enthalpy of Atomization of Reactants ---
    # The formation reaction is: 12 C(s) + 11 H2(g) -> C12H22(g)
    num_H2_moles = num_H / 2
    
    reactant_atomization_energy = (num_C * H_atom_C) + (num_H2_moles * BE_HH)
    # Calculation: (12 * 1000) + (11 * 100) = 12000 + 1100 = 13100 kJ

    # --- 3. Calculate Enthalpy of Atomization of Product (Sum of Bond Energies) ---
    product_bond_energy = (num_CC_single_bonds * BE_CC_single) + \
                          (num_CC_double_bonds * BE_CC_double) + \
                          (num_CH_bonds * BE_CH)
    # Calculation: (9 * 200) + (2 * 300) + (22 * 400) = 1800 + 600 + 8800 = 11200 kJ/mol

    # --- 4. Calculate Enthalpy of Formation (ΔHf) in kJ/mol ---
    delta_Hf_mol = reactant_atomization_energy - product_bond_energy
    # Calculation: 13100 - 11200 = 1900 kJ/mol

    # --- 5. Convert to kJ/g ---
    # Using integer atomic masses (C=12, H=1) as implied by the solutions
    molar_mass = (num_C * 12.0) + (num_H * 1.0)
    # Calculation: (12 * 12) + (22 * 1) = 144 + 22 = 166 g/mol

    delta_Hf_g = delta_Hf_mol / molar_mass
    # Calculation: 1900 / 166 ≈ 11.4457 kJ/g

    # --- 6. Compare with Options and Verify the Answer ---
    # The provided answer is 'A', which corresponds to the value 11.44 kJ/g.
    # The options given in the prompt are:
    # A) 11.44 kJ/g
    # B) 1900 kJ/g
    # C) 11200 kJ/mol
    # D) 67.4 kJ/mol
    
    final_answer_value = 11.44

    # Check if the calculated value matches the value from option A
    if not math.isclose(delta_Hf_g, final_answer_value, rel_tol=1e-2):
        return f"Incorrect. The calculated enthalpy of formation is {delta_Hf_g:.2f} kJ/g, which does not match the value in option A ({final_answer_value} kJ/g)."

    # Verify that other options are indeed distractors
    # Check Option C: 11200 kJ/mol
    if not math.isclose(product_bond_energy, 11200):
        return "Incorrect. The logic for identifying distractors is flawed. The calculated product bond energy is not 11200 kJ/mol."
    
    # Check Option B: 1900 kJ/g
    if not math.isclose(delta_Hf_mol, 1900):
        return "Incorrect. The logic for identifying distractors is flawed. The calculated enthalpy of formation is not 1900 kJ/mol."
    
    # If the calculated value matches option A and the distractors are correctly identified, the answer is correct.
    return "Correct"

# Run the check
result = check_enthalpy_of_formation()
print(result)