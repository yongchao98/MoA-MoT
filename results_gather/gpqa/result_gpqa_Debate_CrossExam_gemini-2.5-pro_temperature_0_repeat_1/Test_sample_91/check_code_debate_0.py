import math

def check_enthalpy_of_formation_answer():
    """
    This function checks the correctness of the calculated enthalpy of formation for
    (CH3)2C=CH-CH2-CH(CH3)-CH2-CH=C(CH3)2.
    It recalculates the value based on the provided data and compares it with the given answer.
    """
    # --- Define problem constraints and given data ---
    enthalpy_atomization_C = 1000  # kJ/mol
    bond_energy_H_H = 100         # kJ/mol
    bond_energy_C_C = 200         # kJ/mol
    bond_energy_C_eq_C = 300      # kJ/mol
    bond_energy_C_H = 400         # kJ/mol
    atomic_mass_C = 12.0          # g/mol
    atomic_mass_H = 1.0           # g/mol

    # --- Step 1: Analyze the molecule's structure ---
    # Molecule: (CH3)2C=CH-CH2-CH(CH3)-CH2-CH=C(CH3)2
    # The interpretation of the chemical formula leads to the following counts,
    # which are correctly identified in the provided solution.
    num_C_atoms = 12
    num_H_atoms = 22
    num_C_H_bonds = 22
    num_C_C_bonds = 9
    num_C_eq_C_bonds = 2

    # --- Step 2: Calculate total enthalpy of atomization of elements ---
    # The compound C12H22 is formed from 12 C(s) and 11 H2(g).
    # Energy to atomize 12 moles of Carbon
    atomization_energy_C = num_C_atoms * enthalpy_atomization_C
    
    # Energy to atomize 11 moles of H2 (breaking 11 moles of H-H bonds)
    num_H2_moles = num_H_atoms / 2
    atomization_energy_H = num_H2_moles * bond_energy_H_H
    
    total_atomization_energy_elements = atomization_energy_C + atomization_energy_H

    # --- Step 3: Calculate total bond energy of the molecule ---
    # This is the energy released when 1 mole of the gaseous molecule is formed from its gaseous atoms.
    bond_energy_from_CH = num_C_H_bonds * bond_energy_C_H
    bond_energy_from_CC = num_C_C_bonds * bond_energy_C_C
    bond_energy_from_CeqC = num_C_eq_C_bonds * bond_energy_C_eq_C
    
    total_bond_energy_molecule = bond_energy_from_CH + bond_energy_from_CC + bond_energy_from_CeqC

    # --- Step 4: Calculate the enthalpy of formation (ΔHf°) in kJ/mol ---
    # ΔHf° = Σ(ΔH_atom of reactants) - Σ(Bond energies of product)
    delta_Hf_mol = total_atomization_energy_elements - total_bond_energy_molecule

    # --- Step 5: Convert to kJ/g ---
    molar_mass = (num_C_atoms * atomic_mass_C) + (num_H_atoms * atomic_mass_H)
    if molar_mass == 0:
        return "Incorrect: Molar mass cannot be zero."
        
    delta_Hf_g = delta_Hf_mol / molar_mass

    # --- Verification against the provided solution ---
    # The provided answer selects option 'C', which is 11.44 kJ/g.
    # Let's check if our independent calculation matches this result.

    # Check intermediate values from the provided solution for a detailed check
    solution_total_atomization_energy_elements = 13100
    solution_total_bond_energy_molecule = 11200
    solution_delta_Hf_mol = 1900
    solution_molar_mass = 166
    
    if not math.isclose(total_atomization_energy_elements, solution_total_atomization_energy_elements):
        return f"Incorrect: The calculation for the total enthalpy of atomization of elements is wrong. The solution states {solution_total_atomization_energy_elements}, but the code calculated {total_atomization_energy_elements}."
        
    if not math.isclose(total_bond_energy_molecule, solution_total_bond_energy_molecule):
        return f"Incorrect: The calculation for the total bond energy of the molecule is wrong. The solution states {solution_total_bond_energy_molecule}, but the code calculated {total_bond_energy_molecule}."

    if not math.isclose(delta_Hf_mol, solution_delta_Hf_mol):
        return f"Incorrect: The calculation for the enthalpy of formation in kJ/mol is wrong. The solution states {solution_delta_Hf_mol}, but the code calculated {delta_Hf_mol}."

    if not math.isclose(molar_mass, solution_molar_mass):
        return f"Incorrect: The calculation for the molar mass is wrong. The solution states {solution_molar_mass}, but the code calculated {molar_mass}."

    # Check the final answer value against option C
    option_C_value = 11.44
    # We use a relative tolerance to account for potential rounding in the option value.
    if not math.isclose(delta_Hf_g, option_C_value, rel_tol=1e-3):
        return f"Incorrect: The final calculated value is {delta_Hf_g:.4f} kJ/g, which does not match option C ({option_C_value} kJ/g). The selected option is wrong."
        
    # If all checks pass, the reasoning and the final answer are correct.
    return "Correct"

# Execute the check and print the result.
result = check_enthalpy_of_formation_answer()
print(result)