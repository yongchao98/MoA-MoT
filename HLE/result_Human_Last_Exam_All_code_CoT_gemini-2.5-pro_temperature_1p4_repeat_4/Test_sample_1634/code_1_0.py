def solve_irreducible_space():
    """
    This function demonstrates why n=2 is the smallest integer for which a reducible
    n-point topological space exists.
    """
    print("Finding the smallest non-negative integer n such that an n-point topological space can be reducible.")
    print("-" * 80)

    print("Step 1: Analyze n=0.")
    print("A 0-point space X = {} has no proper subsets. Thus, it cannot be a union of proper closed subsets. It is irreducible.")
    print("\nStep 2: Analyze n=1.")
    print("A 1-point space X = {p} has only one proper subset: the empty set {}. The empty set is always closed.")
    print("The union of any number of empty sets is still the empty set, not X. So, any 1-point space is irreducible.")
    print("\nStep 3: Analyze n=2.")
    print("Let the space be X = {1, 2}. We will check if a topology exists that makes it reducible.")
    print("A space is reducible if X = Z1 U Z2, where Z1 and Z2 are proper closed subsets.")
    print("To form X as a union, we can try Z1 = {1} and Z2 = {2}.")
    print("This requires a topology where both {1} and {2} are closed sets.")
    print("Consider the discrete topology on X, where every subset is closed.")
    
    # Define the sets for our demonstration
    X = {1, 2}
    Z1 = {1}
    Z2 = {2}

    print("\nWith the discrete topology on X = {1, 2}:")
    print(f"- Z1 = {Z1} is a proper subset of X and is closed.")
    print(f"- Z2 = {Z2} is a proper subset of X and is closed.")
    
    # Calculate the union
    union_set = Z1.union(Z2)

    # Output the equation with each number as requested
    z1_list = sorted(list(Z1))
    z2_list = sorted(list(Z2))
    union_list = sorted(list(union_set))

    print("\nLet's check their union, which forms the equation:")
    print(f"{{{z1_list[0]}}} U {{{z2_list[0]}}} = {{{', '.join(map(str, union_list))}}}")

    print(f"\nThe union {union_set} is equal to the original space X = {X}.")
    print("Since X is a union of two of its proper closed subsets, this 2-point space is reducible.")
    
    print("-" * 80)
    print("Conclusion: The smallest non-negative integer n for which a reducible n-point space exists is 2.")

solve_irreducible_space()