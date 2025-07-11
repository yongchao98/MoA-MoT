import itertools

def main():
    """
    Finds the smallest n for which an n-point space can be non-irreducible.
    """
    print("An irreducible space cannot be written as a finite union of its proper closed subsets.")
    print("A space is not irreducible (or reducible) if X = Z₁ ∪ Z₂ ∪ ... ∪ Zₖ, where Zᵢ are proper closed subsets.")
    print("\nWe are looking for the smallest number of points, n, in a space that can be reducible.")

    # We test small values of n.
    # n=0: The empty set is irreducible.
    # n=1: A one-point space is always irreducible.
    print("\nChecking n = 0 and n = 1: Any space with 0 or 1 points is always irreducible.")
    
    # Check n=2
    n = 2
    print(f"\nChecking n = {n}:")
    points = set(range(n)) # X = {0, 1}

    # Define the discrete topology for X = {0, 1}
    # Open sets are all subsets of X.
    # This means closed sets are also all subsets of X.
    print(f"Let the space X be the set of points: {points}")
    print("Consider the discrete topology, where every subset of X is both open and closed.")

    # We choose two proper closed subsets. In the discrete topology, singleton sets are closed.
    Z1 = {0}
    Z2 = {1}
    
    # Both Z1 and Z2 are proper subsets of X, and they are closed.
    print(f"We can choose two proper closed subsets:")
    print(f"  Z₁ = {Z1}")
    print(f"  Z₂ = {Z2}")
    
    # Calculate the union
    union_of_subsets = Z1.union(Z2)

    # Demonstrate that their union is the whole space X.
    print("\nThe union of these subsets is:")
    
    # Outputting each number in the final equation
    # To format the sets for printing, e.g., "{0, 1}"
    z1_str = "{" + ", ".join(map(str, sorted(list(Z1)))) + "}"
    z2_str = "{" + ", ".join(map(str, sorted(list(Z2)))) + "}"
    union_str = "{" + ", ".join(map(str, sorted(list(union_of_subsets)))) + "}"

    print(f"  {z1_str} ∪ {z2_str} = {union_str}")
    
    if union_of_subsets == points:
        print(f"\nThe union equals the original space X. Therefore, the space is not irreducible.")
        print(f"We have found a non-irreducible space with n = {n} points.")
        print(f"\nThe smallest non-negative integer n is {n}.")
        
        # The final answer in the required format
        # This will be captured by the system.
        print(f"\n<<<{n}>>>")

if __name__ == "__main__":
    main()
