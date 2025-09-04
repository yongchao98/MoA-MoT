import math

def check_enthalpy_of_formation():
    """
    This function checks the correctness of the calculated enthalpy of formation for the given molecule.
    It follows the standard procedure using bond energies and enthalpy of atomization.
    """
    # --- Given Data ---
    H_atom_C = 1000      # Enthalpy of atomization of carbon, kJ/mol
    BE_HH = 100          # Bond energy of H-H, kJ/mol
    BE_CC = 200          # Bond energy of C-C, kJ/mol
    BE_C_dbl_C = 300     # Bond energy of C=C, kJ/mol
    BE_CH = 400          # Bond energy of C-H, kJ/mol

    # --- Step 1: Determine Molecular Formula and Bond Counts ---
    # The molecule is (CH3)2C=CH-CH2-CH(CH3)-CH2-CH=C(CH3)2.
    # Molecular Formula: C12H22
    # Bond Counts:
    # - C-H bonds: 22
    # - C=C bonds: 2
    # - C-C bonds: 9 (Total C-C links are 12-1=11. Subtracting 2 double bonds leaves 9 single bonds)
    num_C = 12
    num_H = 22
    num_CH_bonds = 22
    num_CC_bonds = 9
    num_C_dbl_C_bonds = 2

    # --- Step 2: Calculate Enthalpy of Atomization of Reactants ---
    # The formation reaction is: 12 C(s) + 11 H2(g) -> C12H22(g)
    # Energy needed to form gaseous atoms (12 C(g) + 22 H(g)) from reactants.
    num_H2_moles = num_H / 2
    H_atom_reactants = (num_C * H_atom_C) + (num_H2_moles * BE_HH)
    
    # --- Step 3: Calculate Enthalpy of Atomization of the Product ---
    # This is the sum of all bond energies in the product molecule.
    H_atom_product = (num_CC_bonds * BE_CC) + \
                     (num_C_dbl_C_bonds * BE_C_dbl_C) + \
                     (num_CH_bonds * BE_CH)

    # --- Step 4: Calculate the Enthalpy of Formation (ΔH_f) in kJ/mol ---
    # ΔH_f = (Enthalpy of atomization of reactants) - (Enthalpy of atomization of product)
    delta_H_f_mol = H_atom_reactants - H_atom_product

    # --- Step 5: Convert to kJ/g ---
    # Molar mass of C12H22 (using integer atomic weights as implied by the problem's simple values)
    molar_mass_C = 12
    molar_mass_H = 1
    molar_mass_compound = (num_C * molar_mass_C) + (num_H * molar_mass_H)
    
    delta_H_f_gram = delta_H_f_mol / molar_mass_compound

    # --- Step 6: Verify the Final Answer ---
    # The provided answer is B, which is 11.44 kJ/g.
    expected_value = 11.44

    # Check if the calculated value matches the expected answer within a small tolerance.
    if not math.isclose(delta_H_f_gram, expected_value, rel_tol=1e-2):
        # Check for common errors that match other options
        if math.isclose(delta_H_f_mol, 11200):
             return "Incorrect. The value 11200 kJ/mol (Option D) represents the total bond energy (enthalpy of atomization of the product), not the enthalpy of formation."
        if math.isclose(delta_H_f_gram, 1900):
             return "Incorrect. The value 1900 (Option A) is the correct numerical value for the enthalpy of formation, but the units should be kJ/mol, not kJ/g."
        
        return f"Incorrect. The calculated enthalpy of formation is {delta_H_f_gram:.2f} kJ/g, which does not match the selected answer of {expected_value} kJ/g."

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_enthalpy_of_formation()
print(result)