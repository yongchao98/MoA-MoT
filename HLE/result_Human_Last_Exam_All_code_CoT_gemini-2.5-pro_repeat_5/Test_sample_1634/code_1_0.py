def solve_irreducible_space():
    """
    This function explains and finds the smallest non-negative integer n
    such that there exists an n-point topological space that is not irreducible.
    """
    print("Goal: Find the smallest non-negative integer n for which an n-point topological space exists that is NOT irreducible (i.e., is reducible).")
    print("\nDefinition: A space X is reducible if it can be written as a finite union of proper closed subsets (Z_i != X).")
    print("X = Z_1 U Z_2 U ... U Z_k")
    print("-" * 70)

    print("\nStep 1: Analyzing n = 0")
    print("Let the space be the 0-point space, X = {} (the empty set).")
    print("In any topology on X, the only closed set is X itself.")
    print("A proper closed subset is a closed subset that is not equal to X. There are no such subsets.")
    print("Therefore, the 0-point space cannot be written as a union of proper closed subsets, so it is irreducible.")
    print("-" * 70)

    print("\nStep 2: Analyzing n = 1")
    print("Let the space be a 1-point space, X = {p}.")
    print("The closed sets in any topology on X must be complements of open sets. The only possible closed sets are {} and {p}.")
    print("The only proper closed subset is the empty set, Z = {}.")
    print("Any finite union of the empty set is still the empty set, which is not equal to X.")
    print("Therefore, any 1-point space is irreducible.")
    print("-" * 70)

    print("\nStep 3: Analyzing n = 2")
    print("Let's see if we can construct a 2-point space that is reducible.")
    print("Let the space be X = {'p1', 'p2'}.")
    print("\nLet's define a topology on X. We'll use the 'discrete topology', where every subset is open.")
    print("Open sets: {}, {'p1'}, {'p2'}, {'p1', 'p2'}")
    print("In the discrete topology, every subset is also closed. The closed sets are:")
    print("Closed sets: {}, {'p1'}, {'p2'}, {'p1', 'p2'}")
    print("\nThe proper closed subsets are the closed subsets not equal to X:")
    print("Proper closed subsets: {}, {'p1'}, {'p2'}")
    print("\nNow, we check if X can be written as a union of its proper closed subsets.")
    print("Let's choose Z1 = {'p1'} and Z2 = {'p2'}.")
    print("Both Z1 and Z2 are proper closed subsets of X.")
    
    # Define the sets for the equation
    z1 = "{'p1'}"
    z2 = "{'p2'}"
    x = "{'p1', 'p2'}"
    
    print("\nLet's form the union. The final equation is:")
    print(f"Z1 U Z2 = {z1} U {z2} = {x}")
    print(f"Since the union equals X, the space X is reducible.")
    print("-" * 70)

    print("\nConclusion:")
    print("We have shown that for n=0 and n=1, all spaces are irreducible.")
    print("We have constructed a reducible space for n=2.")
    print("Therefore, the smallest non-negative integer n is 2.")

solve_irreducible_space()
<<<2>>>