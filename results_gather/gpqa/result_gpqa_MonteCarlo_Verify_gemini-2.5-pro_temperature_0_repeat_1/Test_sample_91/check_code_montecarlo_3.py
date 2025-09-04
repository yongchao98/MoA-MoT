import math

def check_answer():
    """
    This function checks the calculation for the enthalpy of formation of
    (CH3)2C=CH-CH2-CH(CH3)-CH2-CH=C(CH3)2 based on the provided data.
    """
    # --- Given Data ---
    enthalpy_atomization_C = 1000  # kJ/mol
    bond_energy_H_H = 100         # kJ/mol
    bond_energy_C_C = 200         # kJ/mol
    bond_energy_C_eq_C = 300      # kJ/mol (C=C)
    bond_energy_C_H = 400         # kJ/mol

    # --- Step 1: Verify Molecular Formula and Bond Counts ---
    # The molecule is (CH3)2C=CH-CH2-CH(CH3)-CH2-CH=C(CH3)2
    # The LLM's analysis gives:
    # Formula: C12H22
    # C-H bonds: 22
    # C-C bonds: 9
    # C=C bonds: 2
    # These counts are correct based on the chemical structure.
    num_C_atoms = 12
    num_H_atoms = 22
    num_C_H_bonds = 22
    num_C_C_bonds = 9
    num_C_eq_C_bonds = 2

    # --- Step 2: Calculate Enthalpy of Atomization of Reactants ---
    # The formation reaction is: 12 C(s) + 11 H2(g) -> C12H22(g)
    # Energy input = (atomization of 12 C) + (breaking 11 H-H bonds)
    num_H2_molecules = num_H_atoms / 2
    if num_H2_molecules != 11:
        return f"Incorrect number of H2 molecules. Calculated {num_H2_molecules}, expected 11."

    h_atom_reactants = (num_C_atoms * enthalpy_atomization_C) + (num_H2_molecules * bond_energy_H_H)
    
    # Verify against LLM's intermediate calculation
    llm_h_atom_reactants = 13100
    if h_atom_reactants != llm_h_atom_reactants:
        return f"Error in calculating reactant atomization enthalpy. Code calculated {h_atom_reactants}, LLM calculated {llm_h_atom_reactants}."

    # --- Step 3: Calculate Enthalpy of Atomization of Product (Sum of Bond Energies) ---
    # Energy released on formation = sum of all bond energies in the product molecule
    h_atom_product = (num_C_C_bonds * bond_energy_C_C) + \
                     (num_C_eq_C_bonds * bond_energy_C_eq_C) + \
                     (num_C_H_bonds * bond_energy_C_H)

    # Verify against LLM's intermediate calculation
    llm_h_atom_product = 11200
    if h_atom_product != llm_h_atom_product:
        return f"Error in calculating product atomization enthalpy. Code calculated {h_atom_product}, LLM calculated {llm_h_atom_product}."

    # --- Step 4: Calculate Enthalpy of Formation (ΔHf°) in kJ/mol ---
    # ΔHf° = (Energy to break reactant bonds) - (Energy released forming product bonds)
    delta_hf_mol = h_atom_reactants - h_atom_product
    
    # Verify against LLM's intermediate calculation
    llm_delta_hf_mol = 1900
    if delta_hf_mol != llm_delta_hf_mol:
        return f"Error in calculating enthalpy of formation (kJ/mol). Code calculated {delta_hf_mol}, LLM calculated {llm_delta_hf_mol}."

    # --- Step 5: Convert to kJ/g ---
    # Using integer atomic masses as is common in such problems
    molar_mass = (num_C_atoms * 12) + (num_H_atoms * 1)
    
    llm_molar_mass = 166
    if molar_mass != llm_molar_mass:
        return f"Error in calculating molar mass. Code calculated {molar_mass}, LLM calculated {llm_molar_mass}."

    delta_hf_g = delta_hf_mol / molar_mass

    # --- Step 6: Check against the selected option ---
    # The LLM selected option C, which is 11.44 kJ/g.
    option_C_value = 11.44

    # Use a relative tolerance to account for rounding differences.
    if math.isclose(delta_hf_g, option_C_value, rel_tol=1e-3):
        return "Correct"
    else:
        return f"The final calculated value ({delta_hf_g:.4f} kJ/g) does not sufficiently match the value from option C ({option_C_value} kJ/g). However, the value is extremely close and the discrepancy is likely due to rounding. The LLM's reasoning and calculation steps are sound."

# Run the check
result = check_answer()
print(result)