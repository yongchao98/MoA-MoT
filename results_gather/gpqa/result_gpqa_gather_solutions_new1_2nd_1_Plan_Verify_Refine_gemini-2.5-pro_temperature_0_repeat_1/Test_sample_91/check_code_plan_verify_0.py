import math

def check_answer():
    """
    This function checks the correctness of the given answer by recalculating the enthalpy of formation.
    """
    # Given thermodynamic data from the question
    enthalpy_atomization_C = 1000  # kJ/mol
    bond_energy_HH = 100          # kJ/mol
    bond_energy_CC = 200          # kJ/mol
    bond_energy_C_dbl_C = 300     # kJ/mol
    bond_energy_CH = 400          # kJ/mol

    # Step 1: Define molecular structure properties
    # Molecule: (CH3)2C=CH-CH2-CH(CH3)-CH2-CH=C(CH3)2
    num_C_atoms = 12
    num_H_atoms = 22
    num_CC_bonds = 9
    num_C_dbl_C_bonds = 2
    num_CH_bonds = 22

    # Step 2: Calculate the enthalpy of atomization of reactants
    # Formation reaction: 12 C(s) + 11 H2(g) -> C12H22(g)
    moles_of_H2 = num_H_atoms / 2
    energy_atomize_reactants = (num_C_atoms * enthalpy_atomization_C) + (moles_of_H2 * bond_energy_HH)

    # Step 3: Calculate the sum of bond energies in the product (enthalpy of atomization of product)
    energy_form_product_bonds = (num_CC_bonds * bond_energy_CC) + \
                                (num_C_dbl_C_bonds * bond_energy_C_dbl_C) + \
                                (num_CH_bonds * bond_energy_CH)

    # Step 4: Calculate the enthalpy of formation in kJ/mol
    # ΔHf = (Energy to atomize reactants) - (Energy released forming product)
    delta_Hf_mol = energy_atomize_reactants - energy_form_product_bonds

    # Step 5: Convert the enthalpy of formation to kJ/g
    # Using integer atomic masses as implied by the problem's simple numbers
    molar_mass = (12 * 12) + (22 * 1)
    delta_Hf_g = delta_Hf_mol / molar_mass

    # The options provided in the question
    options = {
        "A": 67.4,   # kJ/mol
        "B": 1900,   # kJ/g
        "C": 11.44,  # kJ/g
        "D": 11200   # kJ/mol
    }
    
    # The candidate's answer
    candidate_answer_letter = "C"
    candidate_answer_value = options[candidate_answer_letter]

    # Check if the calculated value matches the candidate's answer
    # Using a tolerance for floating point comparison
    if math.isclose(delta_Hf_g, candidate_answer_value, rel_tol=1e-2):
        return "Correct"
    else:
        # Provide a detailed reason for the incorrectness
        reason = f"The candidate's answer is incorrect.\n"
        reason += f"The calculated enthalpy of formation is approximately {delta_Hf_g:.2f} kJ/g.\n"
        reason += f"The candidate chose option {candidate_answer_letter}, which corresponds to {candidate_answer_value} kJ/g.\n"
        reason += "Calculation steps:\n"
        reason += f"1. Energy to atomize reactants (12 C, 11 H2): {energy_atomize_reactants} kJ\n"
        reason += f"2. Sum of bond energies in product (C12H22): {energy_form_product_bonds} kJ/mol\n"
        reason += f"3. Enthalpy of formation (kJ/mol): {energy_atomize_reactants} - {energy_form_product_bonds} = {delta_Hf_mol} kJ/mol\n"
        reason += f"4. Molar mass of C12H22: {molar_mass} g/mol\n"
        reason += f"5. Enthalpy of formation (kJ/g): {delta_Hf_mol} / {molar_mass} = {delta_Hf_g:.2f} kJ/g\n"
        reason += f"The correct option is C."
        return reason

# Run the check
result = check_answer()
print(result)