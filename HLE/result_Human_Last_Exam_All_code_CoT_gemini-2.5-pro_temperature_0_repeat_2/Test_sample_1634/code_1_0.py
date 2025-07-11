def solve_irreducible_space_problem():
    """
    This function demonstrates that a 2-point space can be reducible,
    which helps find the smallest n for a non-irreducible space.
    """
    # For n=2, let the space be X = {1, 2}.
    X = {1, 2}

    # In the discrete topology on X, every subset is closed.
    # We choose two proper subsets of X.
    Z1 = {1}
    Z2 = {2}

    print("Let's analyze the case for n = 2 to find the smallest n for a reducible space.")
    print(f"Let the topological space be X = {X}.")
    print("We consider the discrete topology, where every subset is closed.")
    print(f"Let's choose two proper closed subsets: Z1 = {Z1} and Z2 = {Z2}.")
    
    print("\nA space is reducible if it can be written as the union of its proper closed subsets.")
    print(f"We check if X can be written as Z1 U Z2.")

    # Verify that Z1 and Z2 are proper subsets of X
    is_z1_proper = Z1.issubset(X) and Z1 != X
    is_z2_proper = Z2.issubset(X) and Z2 != X

    # Verify that their union equals X
    union_of_subsets = Z1.union(Z2)
    is_union_equal_to_X = (union_of_subsets == X)

    print(f"\nIs Z1 a proper subset of X? {is_z1_proper}")
    print(f"Is Z2 a proper subset of X? {is_z2_proper}")
    print(f"Is the union Z1 U Z2 equal to X? {is_union_equal_to_X}")

    if is_z1_proper and is_z2_proper and is_union_equal_to_X:
        print("\nConclusion: A 2-point space can be reducible.")
        
        # Format the sets for clear printing of the equation
        z1_str = "{" + ", ".join(map(str, sorted(list(Z1)))) + "}"
        z2_str = "{" + ", ".join(map(str, sorted(list(Z2)))) + "}"
        x_str = "{" + ", ".join(map(str, sorted(list(X)))) + "}"
        
        print(f"The decomposition is shown by the equation: {z1_str} U {z2_str} = {x_str}")
    else:
        # This part of the code should not be reached with the chosen Z1 and Z2
        print("\nThis choice of subsets does not show that the space is reducible.")

    print("\nSince 0-point and 1-point spaces are always irreducible, the smallest non-negative integer n")
    print("such that an n-point space can be reducible is 2.")

solve_irreducible_space_problem()