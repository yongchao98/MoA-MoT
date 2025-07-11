def solve():
    """
    Finds the smallest n for which an n-point space can be reducible.
    The logic has shown that n=0 and n=1 spaces are always irreducible.
    We demonstrate that for n=2, a reducible space exists.
    """

    # For n=2, we consider a set X with two points.
    point1 = 0
    point2 = 1
    X = {point1, point2}

    # With the discrete topology on X, every subset is closed.
    # We choose two subsets Z1 and Z2.
    Z1 = {point1}
    Z2 = {point2}

    # In the discrete topology, Z1 and Z2 are proper closed subsets.
    # A space is reducible if it is a union of its proper closed subsets.
    # We verify that the union of Z1 and Z2 is X.
    is_reducible = (Z1.union(Z2) == X)

    if is_reducible:
        # The prompt requires printing each number in the final equation.
        # Python's set representation e.g. {0, 1} handles this naturally.
        print("For n=2, a space X can be reducible.")
        print("Example: X = {0, 1} with the discrete topology.")
        print("Here, Z1 = {0} and Z2 = {1} are proper closed subsets.")
        print("The equation showing X is reducible is:")
        print(f"{X} = {Z1} U {Z2}")
    
    # The smallest non-negative integer n is 2.
    answer = 2
    print(f"<<<{answer}>>>")

solve()