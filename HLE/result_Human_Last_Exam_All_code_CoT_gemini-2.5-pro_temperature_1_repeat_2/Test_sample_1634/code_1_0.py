def solve_irreducible_space_problem():
    """
    This script finds the smallest n for which an n-point space can be non-irreducible.
    It does this by analyzing n=0, n=1, and n=2.
    """
    
    # Step 1: Explain the core concepts
    print("--- Problem Definition ---")
    print("A topological space X is 'not irreducible' (or 'reducible') if it can be written as a union of a finite number of 'proper' closed subsets.")
    print("A 'proper' subset is any subset that is not the whole space X itself.")
    print("We are looking for the smallest number of points, 'n', for which such a space can be constructed.\n")

    # Step 2: Analyze the smallest cases
    print("--- Analysis of Small Spaces ---")
    print("Case n = 0:")
    print("Let X be the empty set {}. The only closed subset is {} itself. It is not a *proper* subset, so we cannot form a union of proper closed subsets. The 0-point space is irreducible.\n")

    print("Case n = 1:")
    print("Let X = {p1}. The only proper subset of X is the empty set {}. For X to be reducible, it would have to be a union of proper closed subsets. The only candidate for a proper closed subset is {}. Any union of {} is still {}, which is not X. So, any 1-point space is irreducible.\n")

    # Step 3: Analyze n=2 and demonstrate it is the answer
    print("Case n = 2:")
    print("Let's try to construct a non-irreducible space with n=2 points.")
    
    # Define the 2-point set
    X = {'p1', 'p2'}
    print(f"Let our space be X = {X}.")

    print("\nConsider the 'discrete topology' on X. In this topology, every subset is open, which also means every subset is closed.")
    print(f"The closed subsets of X are: {{}}, {{'p1'}}, {{'p2'}}, and {X}.")

    print("\nTo show X is not irreducible, we need to find two *proper* closed subsets that cover X.")
    # Define the two proper closed subsets
    Z1 = {'p1'}
    Z2 = {'p2'}
    
    print(f"Let's choose Z1 = {Z1}. This is a proper closed subset of X.")
    print(f"Let's choose Z2 = {Z2}. This is also a proper closed subset of X.")

    # Step 4: Show the union and print the final equation
    print("\nNow, let's form their union, Z1 U Z2, to see if it equals X.")
    
    # The prompt requires printing each "number" in the final equation.
    # We will represent the elements of the sets.
    p1 = list(Z1)[0]
    p2 = list(Z2)[0]
    union_result = Z1.union(Z2)

    print("\n--- The Final Equation ---")
    print(f"Z1 U Z2 = {{{p1}}} U {{{p2}}}")
    print(f"The result of the union is: {union_result}")
    print(f"This is equal to our original space X = {X}.")
    
    # Step 5: Final Conclusion
    print("\n--- Conclusion ---")
    print("Since we found two proper closed subsets whose union is X, the 2-point space with the discrete topology is NOT irreducible.")
    print("We have shown that spaces with n=0 and n=1 are always irreducible.")
    print("Therefore, the smallest non-negative integer n for which a non-irreducible space exists is 2.")

solve_irreducible_space_problem()