import re
import math

def check_answer():
    """
    This function checks the correctness of the calculated enthalpy of formation.
    """
    # Given values from the question
    enthalpy_atomization_C = 1000  # kJ/mol
    bond_energy_HH = 100          # kJ/mol
    bond_energy_CC = 200          # kJ/mol
    bond_energy_C_eq_C = 300      # kJ/mol
    bond_energy_CH = 400          # kJ/mol

    # Step 1: Determine the molecular formula and bond counts for (CH3)2C=CH-CH2-CH(CH3)-CH2-CH=C(CH3)2
    # Molecular Formula: C12H22
    num_C_atoms = 12
    num_H_atoms = 22
    
    # Bond Counts
    num_CH_bonds = 22
    num_C_eq_C_bonds = 2
    # Total C-C skeleton bonds = num_C_atoms - 1 = 11
    # num_CC_bonds = total skeleton bonds - num_C_eq_C_bonds = 11 - 2 = 9
    num_CC_bonds = 9

    # Step 2: Calculate the enthalpy of atomization of the reactants
    # Formation reaction: 12 C(s) + 11 H2(g) -> C12H22(g)
    # Energy to atomize reactants = (12 * ΔH_atom(C)) + (11 * BE(H-H))
    enthalpy_atomization_reactants = (num_C_atoms * enthalpy_atomization_C) + ((num_H_atoms / 2) * bond_energy_HH)

    # Step 3: Calculate the total bond energy of the product (enthalpy of atomization of the compound)
    # Energy released forming product = (num_CC * BE(C-C)) + (num_C=C * BE(C=C)) + (num_CH * BE(C-H))
    total_bond_energy_product = (num_CC_bonds * bond_energy_CC) + \
                                (num_C_eq_C_bonds * bond_energy_C_eq_C) + \
                                (num_CH_bonds * bond_energy_CH)

    # Step 4: Calculate the enthalpy of formation (ΔH_f) in kJ/mol
    # ΔH_f = (Enthalpy of atomization of reactants) - (Total bond energy of product)
    delta_H_f_mol = enthalpy_atomization_reactants - total_bond_energy_product

    # Step 5: Convert to kJ/g
    # Molar Mass of C12H22 = (12 * 12) + (22 * 1) = 166 g/mol
    molar_mass = (num_C_atoms * 12) + (num_H_atoms * 1)
    delta_H_f_g = delta_H_f_mol / molar_mass

    # The final answer provided by the LLM
    llm_answer_text = "<<<C>>>"
    
    # Define the options from the question
    options = {
        "A": {"value": 1900, "unit": "kJ/g"},
        "B": {"value": 67.4, "unit": "kJ/mol"},
        "C": {"value": 11.44, "unit": "kJ/g"},
        "D": {"value": 11200, "unit": "kJ/mol"}
    }

    # Extract the letter from the LLM's answer
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return f"Invalid answer format. Expected '<<<A>>>', '<<<B>>>', '<<<C>>>', or '<<<D>>>', but got {llm_answer_text}"

    chosen_option_letter = match.group(1)
    chosen_option = options[chosen_option_letter]
    
    # Check if the chosen option's unit matches the calculation
    if chosen_option["unit"] == "kJ/g":
        calculated_value = delta_H_f_g
    elif chosen_option["unit"] == "kJ/mol":
        calculated_value = delta_H_f_mol
    else:
        return f"Unknown unit in chosen option: {chosen_option['unit']}"

    # Compare the calculated value with the chosen option's value using a tolerance
    if not math.isclose(calculated_value, chosen_option["value"], rel_tol=1e-2):
        return (f"Incorrect. The chosen answer is {chosen_option_letter} ({chosen_option['value']} {chosen_option['unit']}), "
                f"but the calculated value is approximately {delta_H_f_g:.2f} kJ/g (or {delta_H_f_mol} kJ/mol). "
                f"The correct option is C because the calculated enthalpy of formation is 1900 kJ/mol, which is {delta_H_f_g:.2f} kJ/g.")

    # Final check on the calculation steps
    if num_C_atoms != 12 or num_H_atoms != 22:
        return f"Incorrect atom count. Calculated C={num_C_atoms}, H={num_H_atoms}. Expected C=12, H=22."
    if num_CC_bonds != 9 or num_C_eq_C_bonds != 2 or num_CH_bonds != 22:
        return f"Incorrect bond count. Calculated C-C={num_CC_bonds}, C=C={num_C_eq_C_bonds}, C-H={num_CH_bonds}. Expected 9, 2, 22 respectively."
    if not math.isclose(enthalpy_atomization_reactants, 13100):
        return f"Incorrect enthalpy of atomization of reactants. Calculated {enthalpy_atomization_reactants}, expected 13100 kJ."
    if not math.isclose(total_bond_energy_product, 11200):
        return f"Incorrect total bond energy of product. Calculated {total_bond_energy_product}, expected 11200 kJ/mol."
    if not math.isclose(delta_H_f_mol, 1900):
        return f"Incorrect enthalpy of formation in kJ/mol. Calculated {delta_H_f_mol}, expected 1900 kJ/mol."
    if not math.isclose(delta_H_f_g, 11.44, rel_tol=1e-2):
        return f"Incorrect enthalpy of formation in kJ/g. Calculated {delta_H_f_g:.2f}, expected ~11.44 kJ/g."

    return "Correct"

# Run the check
result = check_answer()
print(result)