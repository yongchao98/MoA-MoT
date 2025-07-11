import itertools

def get_powerset(s):
    """Helper function to get the powerset of a set."""
    x = len(s)
    # Using chain to flatten the list of tuples from combinations
    powerset = itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(x + 1))
    return {frozenset(subset) for subset in powerset}

def demonstrate_reducibility_for_n2():
    """
    Demonstrates that a 2-point space can be not irreducible.
    """
    print("We are looking for the smallest non-negative integer n such that an n-point space can be reducible.")
    print("A space is reducible (not irreducible) if it's a union of two proper closed subsets.")
    
    print("\n--- Analysis for n=0 and n=1 ---")
    print("For n=0 or n=1, any space is irreducible because there are not enough proper closed subsets to form the entire space.")
    
    print("\n--- Analysis for n=2 ---")
    n = 2
    # Define the set X with n points
    points = {f'p{i+1}' for i in range(n)}
    X = frozenset(points)
    
    print(f"Let's consider an n={n} point space, X = {set(X)}")

    # We will use the discrete topology, where every subset is open.
    topology = get_powerset(points)
    print("Let's use the discrete topology, where every subset is open.")
    
    # In the discrete topology, every subset is also closed.
    print("In this topology, every subset is also closed because the complement of any open set is also open (and thus also closed).")

    closed_sets = topology
    
    # Proper closed subsets are all closed subsets except X itself.
    proper_closed_sets = {s for s in closed_sets if s != X}
    
    print(f"\nThe proper closed subsets are: {[set(s) for s in proper_closed_sets if s]}")

    # We need to find two proper closed subsets Z1, Z2 such that Z1 U Z2 = X
    Z1 = frozenset({'p1'})
    Z2 = frozenset({'p2'})
    
    print(f"\nLet's choose two proper closed subsets: Z1 = {set(Z1)} and Z2 = {set(Z2)}")
    
    # Calculate the union
    union_set = Z1.union(Z2)
    
    # This is the "final equation" showing the result
    print(f"The union is: {set(Z1)} U {set(Z2)} = {set(union_set)}")
    
    if union_set == X:
        print(f"This union equals our original space X. Therefore, this 2-point space is not irreducible.")
        print(f"\nConclusion: The smallest non-negative integer n is 2.")
    else:
        # This part should not be reached with our chosen sets
        print("Could not demonstrate reducibility with the chosen subsets.")

# Run the demonstration
demonstrate_reducibility_for_n2()

# Final Answer
print("\n<<<2>>>")