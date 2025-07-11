def solve_irreducible_space():
    """
    This function determines the smallest non-negative integer n such that there
    exists an n-point topological space that is not irreducible (i.e., is reducible).
    """

    # Explanation of Irreducibility
    print("Definition: A topological space X is irreducible if it cannot be written as a finite union of proper closed subsets.")
    print("A space is reducible (not irreducible) if X = Z_1 U Z_2 U ... U Z_k, where each Z_i is a closed and proper subset of X.")
    print("\nWe are looking for the smallest number of points, n, for which such a reducible space exists.")

    # Case n=0
    print("\n--- Case n = 0 ---")
    print("Let X be the 0-point space (the empty set). The only closed set is X itself.")
    print("Since there are no proper closed subsets, the 0-point space is irreducible.")

    # Case n=1
    print("\n--- Case n = 1 ---")
    print("Let X = {p}, a 1-point space. The only proper closed subset is the empty set {}.")
    print("The union of any number of empty sets is still the empty set, which is not equal to X.")
    print("Therefore, any 1-point space is irreducible.")

    # Case n=2
    print("\n--- Case n = 2 ---")
    print("Let X = {p1, p2}, a 2-point space.")
    print("We check if a topology exists that makes X reducible. Consider the discrete topology, where every subset is open.")
    print("In this topology, every subset is also closed.")
    print("Let's choose two proper closed subsets:")
    print("  - Z1 = {p1}")
    print("  - Z2 = {p2}")
    print("Both Z1 and Z2 are proper subsets of X and are closed in the discrete topology.")
    
    # Representing the sets for the final equation output
    Z1_str = "{p1}"
    Z2_str = "{p2}"
    X_str = "{p1, p2}"
    
    print("\nWe check if their union equals X.")
    print(f"The final equation is: {Z1_str} U {Z2_str} = {X_str}")
    
    print("\nSince X can be expressed as a union of two proper closed subsets, this 2-point space is reducible.")
    
    smallest_n = 2
    print(f"\nConclusion: The smallest non-negative integer n is {smallest_n}.")

solve_irreducible_space()