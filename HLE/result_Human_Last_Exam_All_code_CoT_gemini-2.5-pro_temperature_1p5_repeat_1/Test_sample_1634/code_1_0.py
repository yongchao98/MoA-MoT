def solve_irreducible_space_problem():
    """
    This script finds the smallest non-negative integer n for which an n-point
    topological space exists that is not irreducible.
    
    Definition: A topological space X is not irreducible (or reducible) if it
    can be written as a finite union of proper closed subsets.
    (X = Z1 U Z2 U ... U Zk, where each Zi is closed and Zi is a proper subset of X).
    """

    # Step 1: Analyze n=0 and n=1.
    # For n=0, X is the empty set. The only closed set is X itself. There are no
    # proper closed subsets, so it's irreducible.
    # For n=1, X = {p}. The only proper closed subset is the empty set.
    # Any union of empty sets is still the empty set, which is not X.
    # So, any 1-point space is irreducible.

    # Step 2: Analyze n=2.
    # We will show that a 2-point space can be reducible.
    print("Finding the smallest n for a non-irreducible n-point space.")
    print("----------------------------------------------------------")
    print("Spaces with n=0 and n=1 points are always irreducible.")
    print("Let's test n=2.")
    
    # Consider a 2-point set.
    X = {1, 2}
    print(f"\nConsider the space X = {X} with {len(X)} points.")

    # We use the discrete topology on X. In the discrete topology, every subset is open.
    # Consequently, every subset is also closed (since its complement is open).
    
    # The closed sets are the power set of X: {}, {1}, {2}, {1, 2}.
    # The proper closed subsets of X are those not equal to X itself.
    # Proper closed subsets: {}, {1}, {2}.
    
    # Let's choose two non-empty proper closed subsets.
    Z1 = {1}
    Z2 = {2}
    
    print("We equip X with the discrete topology, where all subsets are closed.")
    print(f"Let's choose two proper closed subsets of X: Z1 = {Z1} and Z2 = {Z2}.")
    
    # Check if their union equals X.
    union_of_subsets = Z1.union(Z2)
    
    print(f"The union of these subsets is Z1 U Z2 = {union_of_subsets}.")
    
    # Conclusion
    if union_of_subsets == X:
        print(f"\nSince X can be expressed as the union of two proper closed subsets,")
        print("the 2-point space with the discrete topology is NOT irreducible.")
        print("\nTherefore, the smallest such n is 2.")
        
        # Final equation demonstration as requested.
        # We need to print the numbers in the final equation.
        # Equation: {1, 2} = {1} U {2}
        
        x_str = "{" + ", ".join(map(str, sorted(list(X)))) + "}"
        z1_str = "{" + ", ".join(map(str, sorted(list(Z1)))) + "}"
        z2_str = "{" + ", ".join(map(str, sorted(list(Z2)))) + "}"
        
        print("\nFinal Equation Demonstration:")
        print(f"{x_str} = {z1_str} U {z2_str}")
    else:
        # This part should not be reached with our choice of Z1, Z2.
        print("Could not demonstrate reducibility for n=2.")

solve_irreducible_space_problem()