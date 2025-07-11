def solve_irreducible_space():
    """
    This function determines the smallest n for which an n-point topological space
    that is not irreducible can exist, and demonstrates the case.
    
    A space X is not irreducible (i.e., reducible) if it can be expressed
    as a finite union of proper closed subsets (Z_i != X), such that X = Z1 U Z2 U ...
    """
    
    # We have determined through analysis that n=0 and n=1 point spaces are always irreducible.
    # We now test n=2.
    n = 2
    
    print(f"The smallest non-negative integer n such that an n-point topological space exists that is not irreducible is {n}.")
    print("\nHere is a demonstration using a 2-point space X = {1, 2} with the discrete topology.")
    print("In the discrete topology, every subset is open, which also means every subset is closed.")
    
    # Define the space
    X = {1, 2}
    
    # Define two proper closed subsets
    # In the discrete topology, {1} and {2} are both closed sets.
    # They are 'proper' because neither is equal to the whole space X.
    Z1 = {1}
    Z2 = {2}
    
    print(f"\nThe space is X = {X}.")
    print(f"Let's choose two proper closed subsets: Z1 = {Z1} and Z2 = {Z2}.")
    
    # Verify that their union equals the original space
    union_of_subsets = Z1.union(Z2)
    
    if union_of_subsets == X:
        print(f"\nThe union of Z1 and Z2 is {union_of_subsets}, which is equal to X.")
        print("Since X can be written as a union of its proper closed subsets, this space is not irreducible.")
        
        # Format the sets for the final equation string
        z1_str = "{" + ", ".join(map(str, sorted(list(Z1)))) + "}"
        z2_str = "{" + ", ".join(map(str, sorted(list(Z2)))) + "}"
        x_str = "{" + ", ".join(map(str, sorted(list(X)))) + "}"
        
        print("\nThe final equation demonstrating this is:")
        print(f"{z1_str} U {z2_str} = {x_str}")

# Run the solver
solve_irreducible_space()