import itertools

def main():
    """
    This script demonstrates that a 2-point space can be reducible,
    making n=2 the smallest integer for which this is possible.
    """
    # The smallest nonnegative integer n is 2.
    # n=0: The empty set is irreducible as it has no proper subsets.
    # n=1: A 1-point space {p} is irreducible. Its only proper subset is the empty
    #      set, and the union of empty sets cannot be {p}.
    # n=2: We can construct a reducible space. Let's demonstrate.

    # 1. Define the 2-point space X.
    X = {0, 1}
    
    # 2. Define a topology. We'll use the discrete topology, where every
    #    subset is an open set. The set of open sets is the power set of X.
    power_set = []
    for i in range(len(X) + 1):
        for subset in itertools.combinations(X, i):
            power_set.append(frozenset(subset))
    
    open_sets = set(power_set)

    # 3. Determine the closed sets. A set is closed if its complement in X is open.
    #    In the discrete topology, since all subsets are open, their complements
    #    are also open, meaning all subsets are also closed.
    closed_sets = open_sets
    
    # 4. A space is reducible if it's a union of finitely many PROPER closed subsets.
    #    A proper subset is any subset that is not equal to the original set X.
    proper_closed_subsets = {s for s in closed_sets if s != X}
    
    print("This script checks if a 2-point space can be reducible (not irreducible).")
    print("-" * 70)
    
    print(f"Let the topological space be X = {set(X)}.")
    print("Let's use the discrete topology, where every subset is open and also closed.")
    print("\nThe proper closed subsets (closed subsets not equal to X) are:")
    # We want a user-friendly printout of the sets
    proper_closed_subsets_printable = [set(s) for s in proper_closed_subsets]
    print(proper_closed_subsets_printable)

    # 5. We look for a finite collection of proper closed subsets whose union is X.
    #    Let's choose Z1 = {0} and Z2 = {1}.
    Z1 = frozenset({0})
    Z2 = frozenset({1})
    
    print(f"\nLet's test if X can be formed by the union of two proper closed subsets: Z1 = {set(Z1)} and Z2 = {set(Z2)}.")
    
    # Verify that Z1 and Z2 are in our list of proper closed subsets.
    if Z1 in proper_closed_subsets and Z2 in proper_closed_subsets:
        print("Both Z1 and Z2 are proper closed subsets of X.")
        
        # Calculate the union.
        union_of_subsets = Z1.union(Z2)
        
        # Check if the union equals X.
        if union_of_subsets == X:
            print(f"The union of Z1 and Z2 is {set(union_of_subsets)}.")
            
            print("\nThis union is equal to the original space X.")
            print("Thus, the space X is reducible because it is the union of its proper closed subsets Z1 and Z2.")
            
            # As requested, printing the equation with each number.
            x_str = f"{{{', '.join(map(str, sorted(list(X))))}}}"
            z1_str = f"{{{', '.join(map(str, sorted(list(Z1))))}}}"
            z2_str = f"{{{', '.join(map(str, sorted(list(Z2))))}}}"
            
            print("\nFinal Equation:")
            print(f"{x_str} = {z1_str} U {z2_str}")
        
    print("-" * 70)
    print("Since a 2-point space can be reducible, and 0-point and 1-point spaces cannot, the smallest n is 2.")

if __name__ == "__main__":
    main()
