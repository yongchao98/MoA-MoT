def calculate_basis_functions():
    """
    Calculates the number of contracted Gaussian functions for toluene (C7H8)
    with a 6-31G basis set.
    """

    # Molecular formula of toluene is C7H8
    num_carbon = 7
    num_hydrogen = 8

    # Number of basis functions per atom for the 6-31G basis set
    # For Carbon (C): 1 core function (1s) + 4 inner valence (2s, 2p) + 4 outer valence (2s, 2p) = 9
    functions_per_carbon = 9
    # For Hydrogen (H): 1 inner valence (1s) + 1 outer valence (1s) = 2
    functions_per_hydrogen = 2

    # Calculate the total number of basis functions
    total_carbon_functions = num_carbon * functions_per_carbon
    total_hydrogen_functions = num_hydrogen * functions_per_hydrogen
    total_functions = total_carbon_functions + total_hydrogen_functions

    # Print the explanation and the calculation
    print("For a 6-31G basis set calculation of toluene (C7H8):")
    print(f"1. Each Carbon (C) atom has {functions_per_carbon} contracted basis functions.")
    print(f"2. Each Hydrogen (H) atom has {functions_per_hydrogen} contracted basis functions.")
    print("\nThe molecule has 7 Carbon atoms and 8 Hydrogen atoms.")
    print("\nThe calculation is:")
    print(f"({num_carbon} C atoms * {functions_per_carbon} functions/C) + ({num_hydrogen} H atoms * {functions_per_hydrogen} functions/H)")
    print(f"= {total_carbon_functions} + {total_hydrogen_functions}")
    print(f"= {total_functions}")

    print("\nThus, the total number of contracted Gaussian functions is:")
    print(total_functions)

if __name__ == '__main__':
    calculate_basis_functions()
<<<79>>>