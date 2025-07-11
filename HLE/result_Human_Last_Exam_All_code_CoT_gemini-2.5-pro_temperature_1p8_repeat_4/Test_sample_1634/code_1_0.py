def find_smallest_non_irreducible_space_size():
    """
    This program demonstrates that the smallest integer n for which an n-point 
    topological space can be reducible (not irreducible) is 2.
    """
    print("This program finds the smallest non-negative integer 'n' for which an n-point topological space exists that is not irreducible.")
    print("An irreducible space cannot be written as a finite union of its proper closed subsets.")
    print("A space that is *not* irreducible is called reducible.\n")

    print("Analysis Summary:")
    print("For n=0: The empty space has no proper closed subsets, so it's irreducible.")
    print("For n=1: A one-point space's only proper closed subset is the empty set. Any union of such sets cannot form the whole space. So, it's irreducible.")
    print("For n=2: We will now construct a reducible 2-point space.\n")

    # Let's model the 2-point space and the discrete topology
    n = 2
    X = {0, 1}
    Z1 = {0}
    Z2 = {1}
    
    print(f"Let the space be X = {X}.")
    print("Consider the discrete topology on X, where every subset is also a closed set.")
    print(f"In this topology, the subset Z1 = {Z1} is a proper closed subset (since it's not equal to X).")
    print(f"The subset Z2 = {Z2} is also a proper closed subset.")
    
    # Perform the union
    union_of_subsets = Z1.union(Z2)
    
    # The user requested to "output each number in the final equation"
    # We will format the set union to look like an equation.
    z1_elements = list(Z1)
    z2_elements = list(Z2)
    union_elements = list(union_of_subsets)
    
    print("\nLet's check if the union of Z1 and Z2 equals X:")
    print(f"Equation: {{{z1_elements[0]}}} U {{{z2_elements[0]}}} = {{{', '.join(map(str, sorted(union_elements)))}}}")
    
    if union_of_subsets == X:
        print(f"\nThe union result {union_of_subsets} is equal to the space X.")
        print("Therefore, X is a union of its proper closed subsets, making it reducible (not irreducible).")
    
    print("\nConclusion: A 2-point space that is not irreducible exists, and since 0- and 1-point spaces are always irreducible, the smallest such non-negative integer n is 2.")

# Execute the demonstration
find_smallest_non_irreducible_space_size()