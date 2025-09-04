import math

def check_answer():
    """
    This function checks the correctness of the final answer for the given chemistry problem.
    It recalculates the enthalpy of formation based on the provided data and compares it
    to the value in the selected option.
    """
    # --- Given Data ---
    H_atom_C = 1000      # Enthalpy of atomization of carbon (kJ/mol)
    BE_HH = 100          # Bond energy of H-H (kJ/mol)
    BE_CC = 200          # Bond energy of C-C (kJ/mol)
    BE_C_dbl_C = 300     # Bond energy of C=C (kJ/mol)
    BE_CH = 400          # Bond energy of C-H (kJ/mol)

    # --- Options from the Question ---
    options = {
        "A": "11200 kJ/mol",
        "B": "1900 kJ/g",
        "C": "67.4 kJ/mol",
        "D": "11.44 kJ/g"
    }
    
    # --- The Final Answer to be Checked ---
    # The provided final answer is <<<D>>>
    final_answer_key = "D"

    # --- Step 1: Determine Molecular Formula and Bond Counts ---
    # The molecule is (CH3)2C=CH-CH2-CH(CH3)-CH2-CH=C(CH3)2.
    # A careful count reveals:
    num_C = 12
    num_H = 22
    # Bond counts:
    num_CH_bonds = 22
    num_CC_bonds = 9
    num_C_dbl_C_bonds = 2

    # --- Step 2: Calculate Enthalpy of Atomization of Reactants ---
    # The formation reaction is: 12 C(s) + 11 H2(g) -> C12H22(g)
    # Energy to form gaseous atoms from reactants: 12 C(g) + 22 H(g)
    num_H2_moles = num_H / 2
    H_atom_reactants = (num_C * H_atom_C) + (num_H2_moles * BE_HH)
    # H_atom_reactants = (12 * 1000) + (11 * 100) = 13100 kJ

    # --- Step 3: Calculate Enthalpy of Atomization of the Product (Sum of Bond Energies) ---
    # This is the energy released when forming the molecule from its gaseous atoms.
    H_atom_product = (num_CH_bonds * BE_CH) + (num_CC_bonds * BE_CC) + (num_C_dbl_C_bonds * BE_C_dbl_C)
    # H_atom_product = (22 * 400) + (9 * 200) + (2 * 300) = 8800 + 1800 + 600 = 11200 kJ/mol

    # --- Step 4: Calculate the Enthalpy of Formation (ΔH_f) in kJ/mol ---
    # ΔH_f = (Enthalpy of atomization of reactants) - (Enthalpy of atomization of product)
    delta_H_f_mol = H_atom_reactants - H_atom_product
    # delta_H_f_mol = 13100 - 11200 = 1900 kJ/mol

    # --- Step 5: Convert Units to kJ/g ---
    # Molar Mass of C12H22 (using integer atomic masses as per the problem's precision)
    molar_mass = (num_C * 12) + (num_H * 1)  # 166 g/mol
    delta_H_f_gram = delta_H_f_mol / molar_mass
    # delta_H_f_gram = 1900 / 166 ≈ 11.4457 kJ/g

    # --- Step 6: Verify the Final Answer ---
    expected_value_str = options[final_answer_key]
    try:
        value_str, unit = expected_value_str.split()
        expected_value = float(value_str)
    except (ValueError, KeyError):
        return f"Error: Could not parse the provided answer option '{final_answer_key}'."

    calculated_value = 0
    if unit == "kJ/mol":
        calculated_value = delta_H_f_mol
    elif unit == "kJ/g":
        calculated_value = delta_H_f_gram
    else:
        return f"Error: Unknown unit '{unit}' in the selected option."

    # Compare the calculated value with the value from the selected option using a tolerance
    if math.isclose(calculated_value, expected_value, rel_tol=1e-2):
        return "Correct"
    else:
        return (f"Incorrect. The final answer is given as {final_answer_key} ({expected_value_str}), but the calculation does not support this.\n"
                f"The calculated enthalpy of formation is {delta_H_f_mol:.2f} kJ/mol, which is equivalent to {delta_H_f_gram:.4f} kJ/g.\n"
                f"The value in option {final_answer_key} is {expected_value} {unit}, which does not match the calculated value {calculated_value:.4f} {unit}.")

# Execute the check and print the result
result = check_answer()
print(result)