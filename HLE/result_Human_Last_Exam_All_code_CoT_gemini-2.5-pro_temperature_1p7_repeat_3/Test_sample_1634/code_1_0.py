def find_smallest_n_for_non_irreducible_space():
    """
    This program demonstrates that the smallest integer n for which an n-point 
    topological space can be not irreducible is 2.
    
    A space is NOT irreducible (or reducible) if it can be written as a finite 
    union of proper closed subsets. This is equivalent to being the union of 
    two proper closed subsets.

    - For n=0 (empty set), there are no proper closed subsets, so it's irreducible.
    - For n=1 (a point), the only proper closed subset is the empty set, so any 
      union of proper closed subsets is empty, not the whole space. It's irreducible.
    - For n=2, we can construct a reducible space. We demonstrate this below.
    """
    
    n = 2
    
    # Define the 2-point set X
    X = frozenset({0, 1})
    
    print(f"We will test if a {n}-point space can be not irreducible.")
    print(f"Let the space be X = {set(X)}.")
    
    # We choose the discrete topology on X. In this topology, every subset is open.
    # Consequently, every subset is also closed, as the complement of any open
    # set is closed.
    print("Let's consider the discrete topology on X, where every subset is closed.")

    # A space is reducible if X = Z1 U Z2, where Z1 and Z2 are proper closed subsets.
    # "Proper" means not equal to X itself.
    
    # In the discrete topology on X = {0, 1}, the subsets {0} and {1} are closed.
    Z1 = frozenset({0})
    Z2 = frozenset({1})
    
    # They are also proper subsets because they are not equal to X.
    print(f"\nConsider the subset Z1 = {set(Z1)}. It is a proper closed subset of X.")
    print(f"Consider the subset Z2 = {set(Z2)}. It is also a proper closed subset of X.")

    # Check if their union equals X.
    union_of_subsets = Z1.union(Z2)

    print(f"\nLet's check if their union is equal to X:")
    print(f"Equation: {set(Z1)} U {set(Z2)} = {set(union_of_subsets)}")
    
    if union_of_subsets == X:
        print("The union is indeed equal to X.")
        print("\nWe have shown that X is a union of two of its proper closed subsets.")
        print("Therefore, this 2-point space is not irreducible.")
        print("\nSince 0-point and 1-point spaces are always irreducible, the smallest")
        print("non-negative integer n for which an n-point space can be")
        print("not irreducible is 2.")
    else:
        # This part of the code will not be reached with the chosen Z1 and Z2.
        print("The construction failed.")

find_smallest_n_for_non_irreducible_space()

<<<2>>>