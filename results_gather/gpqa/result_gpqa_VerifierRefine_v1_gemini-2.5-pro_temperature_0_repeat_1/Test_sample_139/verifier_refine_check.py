import math

def check_correctness():
    """
    This function checks the correctness of the provided LLM answer by verifying each step of the deduction.
    """
    # --- Data Store ---
    # Storing key properties of elements and compounds mentioned in the problem.
    # Z: atomic number, N: neutrons, A: mass number (of the relevant isotope)
    ELEMENT_DATA = {
        'H': {'Z': 1, 'isotopes': {'D': {'N': 1, 'A': 2}}},
        'Li': {'Z': 3, 'isotopes': {'7Li': {'A': 7}}},
        'Al': {'Z': 13, 'isotopes': {'27Al': {'A': 27}}},
    }

    # Storing qualitative chemical properties to verify the narrative.
    CHEMICAL_PROPERTIES = {
        'LiAlD4': {
            'common_analog': 'LiAlH4',
            'is_strong_reducer': True,
            'reacts_violently_with_water': True,
            'formula_map': {'Li': 1, 'Al': 1, 'D': 4}
        },
        'LiAlH4': {
            'is_common_organic_reagent': True
        },
        'Al(OD)3': {
            'is_precipitate': True,
            'decomposition_product_on_heating': 'D2O'
        }
    }

    # --- Verification of the LLM's deductions ---

    # Step 1: Verify Gas W and Substance B
    # Constraint: Gas W's molecule has an equal number of neutrons and protons.
    # Proposed W = D2 (Deuterium gas). A Deuterium atom (D) is a 2H isotope.
    protons_in_D = ELEMENT_DATA['H']['Z']
    neutrons_in_D = ELEMENT_DATA['H']['isotopes']['D']['N']
    protons_in_D2 = 2 * protons_in_D  # 2
    neutrons_in_D2 = 2 * neutrons_in_D  # 2
    if protons_in_D2 != neutrons_in_D2:
        return f"Incorrect: The proposed gas W (D2) does not satisfy the condition of having an equal number of protons and neutrons. Protons: {protons_in_D2}, Neutrons: {neutrons_in_D2}."

    # Constraint: Substance B's melting point is very close to 277 K.
    # Proposed B = D2O (heavy water).
    melting_point_D2O_celsius = 3.82
    melting_point_D2O_kelvin = melting_point_D2O_celsius + 273.15
    given_melting_point_kelvin = 277
    if not math.isclose(given_melting_point_kelvin, melting_point_D2O_kelvin, abs_tol=1.0):
        return f"Incorrect: The melting point of the proposed substance B (D2O), which is {melting_point_D2O_kelvin:.2f} K, is not 'very close' to the given 277 K."

    # Step 2: Verify the proposed reaction and substances X, Y, G
    # Proposed: X=LiAlD4, Y=D2O, G=Al(OD)3
    X, G, B = 'LiAlD4', 'Al(OD)3', 'D2O'
    
    # Constraint: X reacts violently with liquid Y.
    if not CHEMICAL_PROPERTIES[X]['reacts_violently_with_water']:
        return f"Incorrect: The proposed substance X ({X}) is not known to react violently with water/heavy water."

    # Constraint: A precipitate G forms.
    if not CHEMICAL_PROPERTIES[G]['is_precipitate']:
        return f"Incorrect: The proposed product G ({G}) is not a precipitate."

    # Constraint: Heating G releases B.
    if CHEMICAL_PROPERTIES[G]['decomposition_product_on_heating'] != B:
        return f"Incorrect: Heating the proposed precipitate G ({G}) does not release the proposed substance B ({B})."

    # Step 3: Verify the organic chemistry clues
    # Constraint: Reaction of a keto acid with X gives a product with 2 oxygen atoms.
    # LiAlD4 is a strong reducer that reduces keto acids (>=3 oxygens) to diols (2 oxygens).
    if not CHEMICAL_PROPERTIES[X]['is_strong_reducer']:
        return f"Incorrect: The proposed substance X ({X}) is not a strong enough reducer to fulfill the keto acid reaction clue."

    # Constraint: X's analog is a common reagent.
    analog = CHEMICAL_PROPERTIES[X]['common_analog']
    if not CHEMICAL_PROPERTIES.get(analog, {}).get('is_common_organic_reagent', False):
        return f"Incorrect: The analog of X ({analog}) is not a common reagent in organic chemistry."

    # Step 4: Verify the final calculation
    # Substance X is identified as LiAlD4.
    formula = CHEMICAL_PROPERTIES[X]['formula_map'] # {'Li': 1, 'Al': 1, 'D': 4}
    
    # Identify lightest and heaviest elements by atomic number (Z).
    # Elements are H (Z=1), Li (Z=3), Al (Z=13).
    # Lightest is Hydrogen (present as Deuterium). Heaviest is Aluminium.
    
    # Sum the masses of all atoms of the lightest and heaviest elements.
    # Mass of all Hydrogen atoms (as Deuterium, D, mass number 2)
    mass_of_lightest_atoms = formula['D'] * ELEMENT_DATA['H']['isotopes']['D']['A']
    
    # Mass of all Aluminium atoms (as 27Al, mass number 27)
    mass_of_heaviest_atoms = formula['Al'] * ELEMENT_DATA['Al']['isotopes']['27Al']['A']
    
    calculated_cumulative_mass = mass_of_lightest_atoms + mass_of_heaviest_atoms
    
    # The LLM's answer is 35.
    expected_answer = 35
    if calculated_cumulative_mass != expected_answer:
        return f"Incorrect: The final calculation is wrong. The cumulative mass is {calculated_cumulative_mass} (from {mass_of_lightest_atoms} + {mass_of_heaviest_atoms}), not {expected_answer}."

    # All checks passed.
    return "Correct"

# Run the check
result = check_correctness()
print(result)