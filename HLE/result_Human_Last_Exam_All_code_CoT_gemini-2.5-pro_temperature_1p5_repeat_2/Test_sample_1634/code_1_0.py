def solve_irreducible_space_problem():
    """
    Finds and demonstrates the smallest n for which an n-point space can be reducible.
    """
    # Introduction to the problem and definitions
    print("A topological space X is not irreducible (or reducible) if it can be written")
    print("as the union of two proper closed subsets, Z1 and Z2, such that X = Z1 U Z2.")
    print("We are looking for the smallest non-negative integer n for an n-point space")
    print("that can be not irreducible.\n")

    print("For n=0 or n=1, any topological space is irreducible.")
    print("Let's test n=2.\n")

    # Define the 2-point space
    n = 2
    # The points can be represented by any two distinct elements, like 0 and 1.
    X = {0, 1}
    print(f"Consider an n={n} point space, X = {X}.")

    # Define the topology (implicitly, by defining the closed sets)
    # We choose the discrete topology, where every subset is closed.
    print("Let's equip X with the discrete topology, where every subset is a closed set.")

    # Define the two proper closed subsets
    Z1 = {0}
    Z2 = {1}

    # Verify that Z1 and Z2 are proper closed subsets
    print(f"Let Z1 = {Z1}. This is a proper subset of X, and it is closed.")
    print(f"Let Z2 = {Z2}. This is also a proper subset of X, and it is closed.")

    # Show that their union is X
    union_of_subsets = Z1.union(Z2)

    # Format the sets for the final equation output
    # This ensures we "output each number in the final equation"
    z1_str = "{" + ", ".join(map(str, sorted(list(Z1)))) + "}"
    z2_str = "{" + ", ".join(map(str, sorted(list(Z2)))) + "}"
    union_str = "{" + ", ".join(map(str, sorted(list(union_of_subsets)))) + "}"
    x_str = "{" + ", ".join(map(str, sorted(list(X)))) + "}"

    print(f"\nNow let's check the union of Z1 and Z2:")
    print(f"{z1_str} U {z2_str} = {union_str}")
    print(f"This union is equal to the original space X = {x_str}.")

    print("\nSince X can be expressed as a union of two of its proper closed subsets,")
    print("this 2-point space is not irreducible.")
    
    print("\nTherefore, the smallest such non-negative integer is 2.")

solve_irreducible_space_problem()
