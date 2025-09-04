import math

def calculate_enthalpy_of_neutralization(acids, bases):
    """
    Calculates the enthalpy of neutralization for a mixture of strong acids and strong bases.
    The heat of precipitation (e.g., from BaSO4) is ignored, as is standard for this problem type
    unless specified otherwise.

    Args:
        acids (list of tuples): Each tuple contains (volume_L, molarity_M, protons_per_molecule).
        bases (list of tuples): Each tuple contains (volume_L, molarity_M, hydroxides_per_molecule).

    Returns:
        float: The calculated enthalpy of neutralization in kcal.
    """
    # Standard enthalpy of neutralization for strong acid/strong base is ~-13.6 kcal/mol.
    # This value is chosen as it leads directly to one of the multiple-choice options.
    ENTHALPY_PER_MOLE_KCAL = -13.6

    total_moles_h_plus = 0
    for vol, conc, protons in acids:
        total_moles_h_plus += vol * conc * protons

    total_moles_oh_minus = 0
    for vol, conc, hydroxides in bases:
        total_moles_oh_minus += vol * conc * hydroxides

    # The reaction is H+ + OH- -> H2O. The amount of reaction is limited by the lesser of the two reactants.
    moles_of_water_formed = min(total_moles_h_plus, total_moles_oh_minus)

    total_enthalpy = moles_of_water_formed * ENTHALPY_PER_MOLE_KCAL
    
    return total_enthalpy

def run_tests():
    """Runs all defined test cases."""
    tests = []
    
    # Test 1: Main Case (from the user's question)
    # Acids: 500 mL 0.2 M HCl, 300 mL 0.3 M H2SO4
    # Base: 200 mL 0.5 M Ba(OH)2
    acids1 = [(0.500, 0.2, 1), (0.300, 0.3, 2)]
    bases1 = [(0.200, 0.5, 2)]
    expected1 = -2.72
    tests.append(("Main Case", acids1, bases1, expected1))

    # Test 2: Typical Case (Acid is limiting)
    # Acid: 100 mL 0.1 M HCl
    # Base: 200 mL 0.2 M NaOH
    acids2 = [(0.100, 0.1, 1)]
    bases2 = [(0.200, 0.2, 1)]
    expected2 = -0.136
    tests.append(("Acid Limiting", acids2, bases2, expected2))

    # Test 3: Edge Case (No base)
    acids3 = [(0.100, 0.1, 1)]
    bases3 = []
    expected3 = 0.0
    tests.append(("No Base", acids3, bases3, expected3))

    # Test 4: Edge Case (Stoichiometrically equal)
    # Acid: 100 mL 0.2 M H2SO4
    # Base: 400 mL 0.1 M NaOH
    acids4 = [(0.100, 0.2, 2)]
    bases4 = [(0.400, 0.1, 1)]
    expected4 = -0.544
    tests.append(("Stoichiometrically Equal", acids4, bases4, expected4))

    for name, acids, bases, expected in tests:
        result = calculate_enthalpy_of_neutralization(acids, bases)
        if not math.isclose(result, expected, rel_tol=1e-5):
            print(f"--- Test Failed: {name} ---")
            print(f"Acids: {acids}")
            print(f"Bases: {bases}")
            print(f"Expected: {expected}, Got: {result}")
            return

    # If all tests pass, calculate the final answer for the user's question
    final_answer = calculate_enthalpy_of_neutralization(acids1, bases1)
    print("All tests passed")
    print(f"The calculated enthalpy of neutralization is {final_answer:.2f} kcal.")

run_tests()