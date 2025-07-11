def solve_irreducible_space_problem():
    """
    This function demonstrates that a 2-point topological space can be reducible (not irreducible),
    thereby showing that n=2 is the smallest such integer.
    """
    print("--- Analysis of Irreducible Spaces ---")
    
    # A topological space is reducible (not irreducible) if it can be written
    # as a finite union of proper closed subsets. We are looking for the smallest
    # integer 'n' for which an n-point space can be reducible.

    # Case n=0: The empty space X = ∅. It has no proper subsets, so it cannot be a union
    # of proper closed subsets. It is irreducible.
    
    # Case n=1: A one-point space X = {p}. The only proper subset is ∅. In any topology,
    # ∅ is the only proper closed set. The union of ∅ with itself is ∅, not X.
    # So, any 1-point space is irreducible.

    print("Analysis shows that for n=0 and n=1, any topological space is irreducible.")
    print("Now, let's examine the case for n=2.\n")

    # --- Setup for n=2 ---
    n = 2
    # Let the space be X = {0, 1}
    X = {0, 1}
    print(f"Consider an n-point space where n = {n}, so X = {X}.")

    # Define the topology. We'll use the discrete topology.
    # In the discrete topology, every subset is an open set.
    # T = {∅, {0}, {1}, {0, 1}}
    print("Let's give X the discrete topology, where every subset is an open set.")

    # Determine the closed sets. A set is closed if its complement is open.
    # Since every set is open, every set is also closed.
    # Closed sets = {X \ ∅, X \ {0}, X \ {1}, X \ {0, 1}} = {{0, 1}, {1}, {0}, ∅}
    print("In this topology, a set is closed if its complement is open. This means every subset is also a closed set.")
    
    # A proper closed subset is a closed subset that is not equal to X.
    # Proper closed subsets = {{1}, {0}, ∅}
    print(f"The proper closed subsets of X are all closed subsets except for X itself.")

    # Find two proper closed subsets whose union is X.
    Z1 = {0}
    Z2 = {1}
    print(f"Let's select two proper closed subsets: Z1 = {Z1} and Z2 = {Z2}.")
    
    union_Z1_Z2 = Z1.union(Z2)

    # Print the equation and the numbers involved.
    print("\nWe check if their union equals X:")
    print(f"The equation is X = Z1 U Z2")
    z1_num = list(Z1)[0]
    z2_num = list(Z2)[0]
    print(f"Substituting the sets: {{{z1_num}}} U {{{z2_num}}} = {union_Z1_Z2}")

    if union_Z1_Z2 == X:
        print("\nThe union is equal to X. Since X is a union of two of its proper closed subsets, it is NOT irreducible.")
        print("\n--- Conclusion ---")
        print(f"We have found a {n}-point space that is not irreducible. Since spaces with n<2 are always irreducible,")
        print(f"the smallest non-negative integer n is 2.")
    else:
        # This part of the code will not be reached for this example.
        print("The union is not equal to X.")

solve_irreducible_space_problem()