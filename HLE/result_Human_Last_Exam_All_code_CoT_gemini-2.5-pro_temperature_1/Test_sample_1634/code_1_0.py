def solve_irreducible_space():
    """
    This function demonstrates that for n=2, there exists a topological space
    that is not irreducible.
    """
    
    # The smallest non-negative integer is n=2.
    n = 2
    
    # Consider the n-point space X = {1, 2}.
    # We equip it with the discrete topology, where every subset is open.
    # Consequently, every subset is also closed.
    X = {1, 2}

    # In this topology, Z1 = {1} is a proper closed subset.
    Z1 = {1}

    # Z2 = {2} is also a proper closed subset.
    Z2 = {2}

    # The union of these two proper closed subsets is the entire space X.
    union_of_subsets = Z1.union(Z2)

    # We print the result to demonstrate that X is not irreducible.
    print(f"The smallest non-negative integer n for which an n-point space can be not irreducible is {n}.")
    print("\nThis can be shown with a 2-point space X = {1, 2} and the discrete topology.")
    print("In this space, we can find two proper closed subsets whose union is X.")
    
    # Formatting the sets for the final equation output
    z1_str = ", ".join(map(str, sorted(list(Z1))))
    z2_str = ", ".join(map(str, sorted(list(Z2))))
    union_str = ", ".join(map(str, sorted(list(union_of_subsets))))
    
    print("\nLet Z1 = {" + z1_str + "} and Z2 = {" + z2_str + "}.")
    print("The union is Z1 U Z2, which gives the equation:")
    print(f"{{{z1_str}}} U {{{z2_str}}} = {{{union_str}}}")

    print("\nSince X is the union of its proper closed subsets Z1 and Z2, it is not irreducible.")


solve_irreducible_space()
<<<2>>>