def find_smallest_n_for_non_irreducible_space():
    """
    This function determines the smallest integer n for which an n-point 
    non-irreducible topological space exists and demonstrates the proof.
    """
    print("Goal: Find the smallest non-negative integer n for which an n-point space can be non-irreducible.")
    print("Definition: A space X is non-irreducible if X = Z1 U Z2, where Z1 and Z2 are proper closed subsets of X.")
    
    print("\n--- Analysis of small n ---")
    print("Case n=0: X is the empty set. Its only subset is {} itself. There are no proper subsets, so no proper closed subsets. Therefore, it must be irreducible.")
    print("Case n=1: X = {p}. The only proper subset is the empty set {}. If {} is closed, its union is still {}, not X. Therefore, any 1-point space must be irreducible.")
    
    print("\nCase n=2: Let's see if we can construct a non-irreducible space.")
    
    # Define the 2-point space
    n = 2
    X = {1, 2}
    print(f"Let's test n = {n} with the space X = {X}.")

    # To make X reducible, we need to write it as a union of proper closed sets.
    # The only way to do this for X={1,2} is with the sets {1} and {2}.
    # So, we need a topology where both {1} and {2} are closed sets.
    # A set is closed if its complement is open.
    # If {1} is closed, its complement {2} must be open.
    # If {2} is closed, its complement {1} must be open.
    # So, the topology must contain {}, {1}, {2}, and {1, 2} (the discrete topology).
    
    # Define the topology and the sets
    Z1 = {1}
    Z2 = {2}
    
    # In the discrete topology, every subset is open, and thus every subset is also closed.
    print(f"Consider the discrete topology on X, where all subsets are open.")
    print(f"In this topology, all subsets are also closed.")
    
    print(f"\nLet's choose two proper subsets: Z1 = {Z1} and Z2 = {Z2}.")
    
    # Verify they are proper closed subsets
    print(f"Is Z1 a proper closed subset? Yes, it is closed and not equal to X.")
    print(f"Is Z2 a proper closed subset? Yes, it is closed and not equal to X.")
    
    # Verify their union is X
    union_of_sets = Z1.union(Z2)
    print(f"\nNow, let's check their union: {Z1} U {Z2} = {union_of_sets}")
    
    # Final check and conclusion
    if union_of_sets == X:
        print(f"The union is equal to X. We have successfully shown that this {n}-point space is not irreducible.")
        print(f"\nSince n=0 and n=1 spaces are always irreducible, the smallest non-negative integer is 2.")
        
        # Output the equation with each number, as requested.
        print("\nThe equation demonstrating this property is:")
        list_z1 = list(Z1)
        list_z2 = list(Z2)
        list_x = sorted(list(X))
        print(f"{{{list_z1[0]}}} U {{{list_z2[0]}}} = {{{list_x[0]}, {list_x[1]}}}")

find_smallest_n_for_non_irreducible_space()

print("\n<<<2>>>")