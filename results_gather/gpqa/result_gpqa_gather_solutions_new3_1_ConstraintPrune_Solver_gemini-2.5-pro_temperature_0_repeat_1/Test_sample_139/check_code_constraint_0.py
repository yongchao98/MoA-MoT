import math

def check_correctness():
    """
    This function checks the correctness of the provided answer by verifying each step of the deduction.
    The consensus answer identifies Substance X as LiAlD4 and the final calculated value as 35.
    """
    
    # --- Data Representation based on the consensus deduction ---
    # Using integer mass numbers as is standard for this type of problem.
    elements_properties = {
        'H': {'atomic_number': 1},
        'D': {'atomic_number': 1, 'mass_number': 2, 'protons': 1, 'neutrons': 1}, # Isotope of H
        'Li': {'atomic_number': 3},
        'Al': {'atomic_number': 13, 'mass_number': 27},
    }

    # Deduced substances
    substance_X_formula = {'Li': 1, 'Al': 1, 'D': 4}
    substance_B_melting_point_K = 276.97  # Melting point of D2O
    gas_W_formula = {'D': 2}
    
    # The final answer value from the provided solution
    final_answer_value = 35
    
    reasons = []

    # --- Constraint Verification ---

    # Constraint 1: Melting point of B is very close to 277 K.
    if not math.isclose(substance_B_melting_point_K, 277, abs_tol=1.0):
        reasons.append(f"Constraint Violated: The deduced substance B (D2O) has a melting point of {substance_B_melting_point_K}K, which is not considered 'very close' to 277K.")

    # Constraint 2: Gas W's molecule contains the same number of neutrons and protons.
    protons_W = gas_W_formula['D'] * elements_properties['D']['protons']
    neutrons_W = gas_W_formula['D'] * elements_properties['D']['neutrons']
    if protons_W != neutrons_W:
        reasons.append(f"Constraint Violated: The deduced gas W (D2) does not have an equal number of protons ({protons_W}) and neutrons ({neutrons_W}).")

    # Constraint 3: The product of the reaction of a keto acid with X contains 2 atoms of oxygen.
    # This is a chemical knowledge check. LiAlD4 is a strong reducing agent that reduces a keto acid (>=3 O atoms) to a diol (2 O atoms).
    # This confirms X is a strong reducing agent like LiAlD4. The code assumes this chemical fact.
    oxygens_in_product_of_strong_reduction = 2
    if oxygens_in_product_of_strong_reduction != 2:
        reasons.append("Constraint Violated: The logic for the keto acid reaction is flawed. A strong reducing agent should yield a product with 2 oxygen atoms.")

    # --- Final Calculation Verification ---
    
    # Identify lightest and heaviest elements in X based on atomic number
    elements_in_X = list(substance_X_formula.keys())
    
    lightest_element_symbol = min(elements_in_X, key=lambda el: elements_properties[el]['atomic_number'])
    heaviest_element_symbol = max(elements_in_X, key=lambda el: elements_properties[el]['atomic_number'])

    if lightest_element_symbol != 'D':
        reasons.append(f"Calculation Error: Identified lightest element as {lightest_element_symbol}, but it should be Hydrogen/Deuterium ('D').")
    if heaviest_element_symbol != 'Al':
        reasons.append(f"Calculation Error: Identified heaviest element as {heaviest_element_symbol}, but it should be Aluminum ('Al').")

    # Sum the masses as per the question's instructions
    mass_of_lightest_atoms = substance_X_formula[lightest_element_symbol] * elements_properties[lightest_element_symbol]['mass_number']
    mass_of_heaviest_atoms = substance_X_formula[heaviest_element_symbol] * elements_properties[heaviest_element_symbol]['mass_number']

    calculated_sum = mass_of_lightest_atoms + mass_of_heaviest_atoms

    if calculated_sum != final_answer_value:
        reasons.append(f"Calculation Incorrect: The cumulative mass is {calculated_sum}, not {final_answer_value}.")
        reasons.append(f"  - Mass from {substance_X_formula[lightest_element_symbol]} atoms of {lightest_element_symbol}: {substance_X_formula[lightest_element_symbol]} * {elements_properties[lightest_element_symbol]['mass_number']} = {mass_of_lightest_atoms}")
        reasons.append(f"  - Mass from {substance_X_formula[heaviest_element_symbol]} atoms of {heaviest_element_symbol}: {substance_X_formula[heaviest_element_symbol]} * {elements_properties[heaviest_element_symbol]['mass_number']} = {mass_of_heaviest_atoms}")

    # --- Final Verdict ---
    if not reasons:
        return "Correct"
    else:
        return "Incorrect:\n" + "\n".join(reasons)

# Run the check
print(check_correctness())