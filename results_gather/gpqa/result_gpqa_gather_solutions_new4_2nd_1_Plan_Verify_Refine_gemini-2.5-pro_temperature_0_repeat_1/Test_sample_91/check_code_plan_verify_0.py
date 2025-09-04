import math

def check_answer():
    """
    Checks the correctness of the calculated enthalpy of formation.
    """
    # --- Given Data from the Question ---
    enthalpy_atomization_C = 1000  # kJ/mol
    bond_energy_H_H = 100         # kJ/mol
    bond_energy_C_C = 200         # kJ/mol
    bond_energy_C_C_double = 300  # kJ/mol
    bond_energy_C_H = 400         # kJ/mol
    
    # Approximate atomic masses for molar mass calculation
    atomic_mass_C = 12  # g/mol
    atomic_mass_H = 1   # g/mol

    # --- Step 1: Determine Molecular Formula and Bond Count ---
    # Molecule: (CH3)2C=CH-CH2-CH(CH3)-CH2-CH=C(CH3)2
    # Carbon atoms: 2+1+1+1+1(in CH3)+1+1+1+2 = 12
    # Hydrogen atoms: 6+1+2+1+3(in CH3)+2+1+6 = 22
    num_C = 12
    num_H = 22
    
    # Bond counts:
    num_C_H_bonds = 22
    num_C_C_double_bonds = 2
    # For an acyclic molecule, total C-C linkages = num_C - 1
    # C-C single bonds = (total C-C linkages) - (C=C double bonds)
    num_C_C_single_bonds = (num_C - 1) - num_C_C_double_bonds

    # --- Step 2: Calculate Enthalpy of Atomization of Reactants ---
    # Formation reaction: 12 C(s) + 11 H2(g) -> C12H22(g)
    # Energy to atomize 12 moles of C(s)
    energy_atomize_C = num_C * enthalpy_atomization_C
    # Energy to break 11 moles of H2(g) (since we need 22 H atoms)
    num_H2_moles = num_H / 2
    energy_atomize_H2 = num_H2_moles * bond_energy_H_H
    total_reactant_atomization_energy = energy_atomize_C + energy_atomize_H2

    # --- Step 3: Calculate Total Bond Energy of the Product ---
    energy_from_C_H = num_C_H_bonds * bond_energy_C_H
    energy_from_C_C_single = num_C_C_single_bonds * bond_energy_C_C
    energy_from_C_C_double = num_C_C_double_bonds * bond_energy_C_C_double
    total_product_bond_energy = energy_from_C_H + energy_from_C_C_single + energy_from_C_C_double

    # --- Step 4: Calculate Enthalpy of Formation (ΔH_f) in kJ/mol ---
    # ΔH_f = (Energy to atomize reactants) - (Energy released forming product)
    delta_H_formation_mol = total_reactant_atomization_energy - total_product_bond_energy

    # --- Step 5: Convert to kJ/g ---
    molar_mass = (num_C * atomic_mass_C) + (num_H * atomic_mass_H)
    if molar_mass == 0:
        return "Error: Molar mass is zero."
    delta_H_formation_g = delta_H_formation_mol / molar_mass

    # --- Step 6: Verify the Answer ---
    # The final answer provided is <<<A>>>, which corresponds to 11.44 kJ/g.
    expected_value = 11.44
    
    # Check if the calculated value is close to the expected value
    # Using a tolerance to account for potential rounding differences
    if math.isclose(delta_H_formation_g, expected_value, rel_tol=1e-2):
        return "Correct"
    else:
        # Provide a detailed reason for the discrepancy
        reason = (
            f"The final answer is incorrect.\n"
            f"The chosen answer is {expected_value} kJ/g.\n"
            f"However, the calculated value is {delta_H_formation_g:.4f} kJ/g.\n"
            f"Calculation breakdown:\n"
            f"  - Reactant atomization energy: {total_reactant_atomization_energy} kJ\n"
            f"  - Product bond energy: {total_product_bond_energy} kJ/mol\n"
            f"  - Enthalpy of formation: {delta_H_formation_mol} kJ/mol\n"
            f"  - Molar mass: {molar_mass} g/mol\n"
            f"  - Enthalpy of formation (kJ/g): {delta_H_formation_mol} / {molar_mass} = {delta_H_formation_g:.4f} kJ/g"
        )
        return reason

# Run the check
result = check_answer()
print(result)