def find_smallest_reducible_space():
    """
    This function demonstrates that a 2-point space can be reducible,
    confirming that 2 is the smallest n for which this is possible.
    """
    n = 2
    # Use frozensets for hashable sets
    X = frozenset(range(n))

    print(f"Checking for n = {n}")
    # Represent the space X as a regular set for printing
    print(f"Space X = {set(X)}")

    # For a 2-point space, we can construct the discrete topology.
    # In the discrete topology, every subset is an open set.
    # The open sets are the power set of X.
    open_sets = {
        frozenset(),
        frozenset({0}),
        frozenset({1}),
        frozenset({0, 1})
    }
    print(f"Using the discrete topology. Open sets: {[set(s) for s in open_sets]}")

    # A set is closed if its complement is open.
    # In the discrete topology, every subset is also closed.
    closed_sets = {X.difference(s) for s in open_sets}
    print(f"Corresponding closed sets: {[set(s) for s in closed_sets]}")

    # A proper subset of X is any subset except X itself.
    proper_closed_subsets = {s for s in closed_sets if s != X}
    print(f"Proper closed subsets: {[set(s) for s in proper_closed_subsets]}")

    # A space is reducible if it's a union of two proper closed subsets.
    # Let's find Z1 and Z2 from the proper_closed_subsets.
    Z1 = frozenset({0})
    Z2 = frozenset({1})

    # Check if they meet the criteria
    if Z1 in proper_closed_subsets and Z2 in proper_closed_subsets:
        if Z1.union(Z2) == X:
            print("\nFound a decomposition of X into two proper closed subsets:")
            # Format the sets for the final equation printout
            x_str = "{" + ", ".join(sorted(map(str, X))) + "}"
            z1_str = "{" + ", ".join(sorted(map(str, Z1))) + "}"
            z2_str = "{" + ", ".join(sorted(map(str, Z2))) + "}"
            print(f"{x_str} = {z1_str} U {z2_str}")
            print("\nThis shows that a 2-point space can be reducible.")
            print("Since 0-point and 1-point spaces are always irreducible,")
            print("the smallest such non-negative integer is 2.")

find_smallest_reducible_space()
<<<2>>>