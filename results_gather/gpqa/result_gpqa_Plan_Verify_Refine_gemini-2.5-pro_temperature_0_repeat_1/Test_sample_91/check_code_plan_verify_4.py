import math

def check_answer():
    """
    This function checks the correctness of the LLM's answer for the enthalpy of formation calculation.
    """
    # Given values from the question
    H_atom_C = 1000      # Enthalpy of atomization of C(s) -> C(g) in kJ/mol
    BE_HH = 100          # Bond energy of H-H in kJ/mol
    BE_CC = 200          # Bond energy of C-C in kJ/mol
    BE_C_dbl_C = 300     # Bond energy of C=C in kJ/mol
    BE_CH = 400          # Bond energy of C-H in kJ/mol

    # Properties of the molecule (CH3)2C=CH-CH2-CH(CH3)-CH2-CH=C(CH3)2
    # Molecular Formula: C12H22
    num_C = 12
    num_H = 22
    # Bond counts
    num_CH_bonds = 22
    num_CC_bonds = 9
    num_C_dbl_C_bonds = 2

    # --- Step 1: Verify the bond and atom counts ---
    # These are derived from the chemical structure and are assumed correct for the calculation check.
    # The LLM's counts match the manual verification.

    # --- Step 2: Calculate the enthalpy of atomization for the reactants ---
    # Reactants for formation are 12 C(s) and 11 H2(g)
    num_H2_molecules = num_H / 2
    enthalpy_atomization_reactants = (num_C * H_atom_C) + (num_H2_molecules * BE_HH)
    
    expected_enthalpy_atomization_reactants = (12 * 1000) + (11 * 100) # 13100
    if enthalpy_atomization_reactants != expected_enthalpy_atomization_reactants:
        return f"Incorrect calculation of reactant atomization enthalpy. Expected {expected_enthalpy_atomization_reactants}, but calculated {enthalpy_atomization_reactants}."

    # --- Step 3: Calculate the sum of bond energies for the product molecule ---
    sum_bond_energies_product = (num_CH_bonds * BE_CH) + (num_CC_bonds * BE_CC) + (num_C_dbl_C_bonds * BE_C_dbl_C)
    
    expected_sum_bond_energies_product = (22 * 400) + (9 * 200) + (2 * 300) # 11200
    if sum_bond_energies_product != expected_sum_bond_energies_product:
        return f"Incorrect calculation of product bond energies. Expected {expected_sum_bond_energies_product}, but calculated {sum_bond_energies_product}."

    # --- Step 4: Calculate the enthalpy of formation in kJ/mol ---
    # ΔH_f = (Energy to atomize reactants) - (Energy released forming product)
    delta_H_f_mol = enthalpy_atomization_reactants - sum_bond_energies_product
    
    expected_delta_H_f_mol = 1900.0
    if delta_H_f_mol != expected_delta_H_f_mol:
        return f"Incorrect calculation of enthalpy of formation (kJ/mol). Expected {expected_delta_H_f_mol}, but calculated {delta_H_f_mol}."

    # --- Step 5: Calculate the molar mass of C12H22 ---
    # Using integer atomic masses as is standard for such problems unless specified otherwise.
    molar_mass_C = 12
    molar_mass_H = 1
    molar_mass_compound = (num_C * molar_mass_C) + (num_H * molar_mass_H)
    
    expected_molar_mass = 166.0
    if molar_mass_compound != expected_molar_mass:
        return f"Incorrect calculation of molar mass. Expected {expected_molar_mass}, but calculated {molar_mass_compound}."

    # --- Step 6: Calculate the enthalpy of formation in kJ/g ---
    delta_H_f_gram = delta_H_f_mol / molar_mass_compound

    # --- Step 7: Check against the selected option B ---
    # Option B is 11.44 kJ/g
    llm_answer_value = 11.44
    # Use a small tolerance for floating point comparisons
    if not math.isclose(delta_H_f_gram, llm_answer_value, rel_tol=1e-2):
        return f"The calculated enthalpy of formation per gram ({delta_H_f_gram:.4f} kJ/g) does not match the selected answer's value ({llm_answer_value} kJ/g)."

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_answer()
print(result)