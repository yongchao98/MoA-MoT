def solve():
    """
    This function demonstrates that a 2-point space can be reducible,
    making n=2 the smallest integer for which a non-irreducible space exists.
    """
    # Let X be a 2-point set.
    X = {'a', 'b'}

    # Consider the discrete topology on X, where every subset is open.
    # The open sets are: {}, {a}, {b}, {a, b}
    # A set is closed if its complement is open.
    # The closed sets are therefore: {a, b}, {b}, {a}, {}

    # A proper closed subset is a closed set that is not the entire space X.
    # The proper closed subsets are: {b}, {a}, {}
    
    # A space is reducible (not irreducible) if it can be written as the
    # union of two proper closed subsets.
    
    # Let's choose two proper closed subsets:
    Z1 = {'a'}
    Z2 = {'b'}

    # We verify their union equals the original space X.
    union_of_subsets = Z1.union(Z2)

    print("We are looking for the smallest n for which an n-point space can be reducible.")
    print("For n=0 and n=1, all spaces are irreducible.")
    print("\nLet's test n=2 with the space X = {'a', 'b'} and the discrete topology.")
    print("In this topology, the sets Z1 = {'a'} and Z2 = {'b'} are both proper closed subsets.")
    print("\nWe check if their union equals X:")
    
    # The final equation, showing each element.
    # To print each "number" (element) in the equation, we format the sets.
    x_str = f"{{'{sorted(list(X))[0]}', '{sorted(list(X))[1]}'}}"
    z1_str = f"{{'{list(Z1)[0]}'}}"
    z2_str = f"{{'{list(Z2)[0]}'}}"
    
    print(f"{x_str} = {z1_str} U {z2_str}")

    print("\nSince the 2-point space X can be written as the union of two proper closed subsets, it is not irreducible.")
    print("Thus, the smallest such n is 2.")

solve()
<<<2>>>