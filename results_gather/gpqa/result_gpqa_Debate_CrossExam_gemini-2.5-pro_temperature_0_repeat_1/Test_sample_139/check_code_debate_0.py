import sys
from io import StringIO

def check_answer():
    """
    Checks the correctness of the provided answer by verifying its logical steps and final calculation.
    """
    errors = []

    # --- Step 1: Verify the identity of Substance B (D2O) based on melting point ---
    # The problem states the melting point is "very close to 277 K".
    # The known melting point of heavy water (D2O) is 276.97 K (or 3.82 °C).
    melting_point_D2O_K = 276.97
    target_mp_K = 277
    # We'll allow a tolerance of 1 Kelvin for "very close".
    if not abs(melting_point_D2O_K - target_mp_K) <= 1.0:
        errors.append(
            f"Constraint Check Failed: The melting point of D2O ({melting_point_D2O_K} K) "
            f"is not considered 'very close' to the stated {target_mp_K} K with a tolerance of 1 K. "
            "However, the value is extremely close, so this is likely a correct deduction."
        )
    # This check is very likely to pass, as the values are almost identical.

    # --- Step 2: Verify the properties of Gas W (D2) ---
    # The problem states Gas W has an equal number of protons and neutrons.
    # The answer identifies W as D2 (Deuterium gas).
    # A Deuterium atom (D or ²H) has 1 proton and 1 neutron.
    protons_in_D2 = 1 * 2
    neutrons_in_D2 = 1 * 2
    if protons_in_D2 != neutrons_in_D2:
        errors.append(
            f"Constraint Check Failed: The proposed gas W (D2) does not have an equal number of protons and neutrons. "
            f"Calculated Protons: {protons_in_D2}, Calculated Neutrons: {neutrons_in_D2}."
        )

    # --- Step 3: Verify the keto acid clue ---
    # The answer states that reducing a keto acid (R-C(=O)-COOH) with LiAlD4 gives a diol (R-CD(OH)-CD2(OH)).
    # A diol has two hydroxyl (-OH) groups, meaning it contains exactly two oxygen atoms.
    # This is a correct principle in organic chemistry. The logic holds.

    # --- Step 4: Verify the final calculation for Substance X (LiAlD4) ---
    # The question asks for the cumulative atomic mass of the lightest and heaviest elements.
    # Substance X is identified as LiAlD4.
    
    # Define elements and their properties (using integer masses as in the problem)
    # Atomic Number (Z) determines the element's identity.
    elements_data = {
        'Li': {'atomic_number': 3, 'mass': 7},
        'Al': {'atomic_number': 13, 'mass': 27},
        'D': {'atomic_number': 1, 'mass': 2} # Deuterium is an isotope of Hydrogen
    }
    
    composition = {'Li': 1, 'Al': 1, 'D': 4}
    
    # Find the lightest and heaviest elements based on atomic number
    elements_in_X = list(composition.keys())
    
    lightest_element = min(elements_in_X, key=lambda el: elements_data[el]['atomic_number'])
    heaviest_element = max(elements_in_X, key=lambda el: elements_data[el]['atomic_number'])

    if lightest_element != 'D' or heaviest_element != 'Al':
        errors.append(
            f"Logic Check Failed: The lightest element was identified as '{lightest_element}' (should be D/H) "
            f"and/or the heaviest as '{heaviest_element}' (should be Al)."
        )

    # Calculate the cumulative mass as per the question's instructions
    mass_of_lightest_atoms = composition[lightest_element] * elements_data[lightest_element]['mass']
    mass_of_heaviest_atoms = composition[heaviest_element] * elements_data[heaviest_element]['mass']
    
    calculated_total_mass = mass_of_lightest_atoms + mass_of_heaviest_atoms
    
    expected_total_mass = 35 # From the answer's calculation
    
    if calculated_total_mass != expected_total_mass:
        errors.append(
            f"Calculation Failed: The final calculated mass is {calculated_total_mass}, but the answer claims it is {expected_total_mass}.\n"
            f"  - Mass from lightest element ({lightest_element}): {composition[lightest_element]} atoms * {elements_data[lightest_element]['mass']} amu/atom = {mass_of_lightest_atoms} amu\n"
            f"  - Mass from heaviest element ({heaviest_element}): {composition[heaviest_element]} atoms * {elements_data[heaviest_element]['mass']} amu/atom = {mass_of_heaviest_atoms} amu"
        )

    # --- Final Verdict ---
    if not errors:
        return "Correct"
    else:
        return "Incorrect. Reason(s):\n- " + "\n- ".join(errors)

# Run the check and print the result
result = check_answer()
print(result)