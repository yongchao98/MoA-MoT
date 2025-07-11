def calculate_basis_functions():
    """
    Calculates the number of contracted Gaussian functions for toluene (C7H8)
    with the 6-31G basis set.
    """
    # Number of atoms in toluene (C7H8)
    num_carbon = 7
    num_hydrogen = 8

    # Number of contracted Gaussian functions per atom for the 6-31G basis set
    # For Carbon: 1 (core) + 4 (inner valence) + 4 (outer valence) = 9
    funcs_per_carbon = 9
    # For Hydrogen: 1 (inner valence) + 1 (outer valence) = 2
    funcs_per_hydrogen = 2

    # Calculate the total number of functions from Carbon atoms
    total_funcs_carbon = num_carbon * funcs_per_carbon

    # Calculate the total number of functions from Hydrogen atoms
    total_funcs_hydrogen = num_hydrogen * funcs_per_hydrogen

    # Calculate the total number of functions for the molecule
    total_funcs = total_funcs_carbon + total_funcs_hydrogen

    # Print the step-by-step calculation
    print(f"Molecule: Toluene (C{num_carbon}H{num_hydrogen})")
    print("Basis Set: 6-31G\n")
    print("Calculation:")
    print(f"Number of Carbon atoms: {num_carbon}")
    print(f"Number of functions per Carbon atom: {funcs_per_carbon}")
    print(f"Contribution from Carbon = {num_carbon} * {funcs_per_carbon} = {total_funcs_carbon}\n")
    print(f"Number of Hydrogen atoms: {num_hydrogen}")
    print(f"Number of functions per Hydrogen atom: {funcs_per_hydrogen}")
    print(f"Contribution from Hydrogen = {num_hydrogen} * {funcs_per_hydrogen} = {total_funcs_hydrogen}\n")
    print("Total number of contracted Gaussian functions:")
    print(f"{total_funcs_carbon} (from C) + {total_funcs_hydrogen} (from H) = {total_funcs}")

calculate_basis_functions()