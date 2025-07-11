def calculate_contracted_functions():
    """
    Calculates the total number of contracted Gaussian functions for toluene (C7H8)
    with a 6-31G basis set.
    """
    # Number of atoms in toluene (C7H8)
    num_c_atoms = 7
    num_h_atoms = 8

    # Number of contracted Gaussian functions per atom for the 6-31G basis set.
    # For Carbon (1s core, 2s/2p valence): 1 (core) + 4*2 (valence) = 9
    funcs_per_c = 9
    # For Hydrogen (1s valence): 1*2 (valence) = 2
    funcs_per_h = 2

    # Calculate the total number of functions from each element
    total_funcs_c = num_c_atoms * funcs_per_c
    total_funcs_h = num_h_atoms * funcs_per_h

    # Calculate the grand total
    total_funcs = total_funcs_c + total_funcs_h

    # Print the step-by-step calculation
    print("To find the total number of contracted Gaussian functions for Toluene (C7H8) with a 6-31G basis set:")
    print(f"1. A Carbon atom has {funcs_per_c} functions.")
    print(f"2. A Hydrogen atom has {funcs_per_h} functions.")
    print("\nFinal Calculation:")
    print(f"({num_c_atoms} C atoms * {funcs_per_c} funcs/C) + ({num_h_atoms} H atoms * {funcs_per_h} funcs/H) = {total_funcs_c} + {total_funcs_h} = {total_funcs}")

calculate_contracted_functions()