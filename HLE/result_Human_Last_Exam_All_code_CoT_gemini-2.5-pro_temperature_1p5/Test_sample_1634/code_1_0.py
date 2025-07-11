def main():
    """
    This script finds the smallest non-negative integer n such that there exists
    an n-point topological space that is not irreducible.
    """
    # Based on the analysis, the smallest integer n is 2.
    # This script will demonstrate why a 2-point space can be not irreducible.
    n = 2

    # A topological space X is not irreducible (or reducible) if it can be written
    # as a finite union of proper closed subsets.
    # Example: X = Z1 U Z2, where Z1 and Z2 are closed, Z1 != X, and Z2 != X.

    # Let the space X be a set of n points.
    X = set(range(n))  # X = {0, 1}

    # Consider the discrete topology on X. In this topology, every subset is
    # both open and closed.

    # We can choose two subsets, Z1 and Z2.
    Z1 = {0}
    Z2 = {1}

    # In the discrete topology on X = {0, 1}:
    # - Z1 is a proper subset of X.
    # - Z1 is a closed set.
    # - Z2 is a proper subset of X.
    # - Z2 is a closed set.
    
    # We check if their union equals the whole space X.
    union_of_subsets = Z1.union(Z2)

    print(f"We are looking for the smallest n where an n-point space can be not irreducible.")
    print("For n=0 and n=1, any space is irreducible.")
    print(f"\nLet's test n = {n}.")
    print(f"Consider the space X = {sorted(list(X))} with the discrete topology.")
    print("In this topology, every subset is a closed set.")
    print("We can choose the following two proper closed subsets:")
    print(f"  Z1 = {sorted(list(Z1))}")
    print(f"  Z2 = {sorted(list(Z2))}")
    print("\nThe union of Z1 and Z2 is the entire space X:")
    
    # The final equation with each number printed explicitly
    print(f"{sorted(list(union_of_subsets))} = {sorted(list(Z1))} U {sorted(list(Z2))}")

    if union_of_subsets == X:
        print("\nSince X can be expressed as a union of its proper closed subsets, it is not irreducible.")
        print(f"Therefore, the smallest n is {n}.")
    
if __name__ == "__main__":
    main()

<<<2>>>