import collections

def solve():
    """
    Finds and demonstrates the smallest non-negative integer n such that
    an n-point topological space can be not irreducible.
    """

    print("--- Analysis of Irreducibility based on Number of Points (n) ---")

    # Case n=0: The empty set is irreducible as it has no proper subsets.
    print("\nCase n=0: The empty space.")
    print("A space is 'not irreducible' if X = Z1 U Z2 U ..., where Zi are proper closed subsets.")
    print("The empty space has no proper subsets, so it cannot be such a union. It is irreducible.")

    # Case n=1: A one-point space is irreducible.
    print("\nCase n=1: A space X = {p}.")
    print("The only proper subset is the empty set, which is always closed.")
    print("The union of any number of empty sets is still the empty set, not X.")
    print("Therefore, any one-point space is irreducible.")

    # Case n=2: We can construct a space that is not irreducible.
    print("\nCase n=2: A space X = {a, b} with the discrete topology.")
    print("In the discrete topology, every subset is open, so every subset is also closed.")
    
    # Define the 2-point space X
    X = {'a', 'b'}

    # Define two proper closed subsets.
    # In the discrete topology, {a} and {b} are closed sets.
    # They are 'proper' because neither is equal to the full space X.
    Z1 = {'a'}
    Z2 = {'b'}

    # The union of Z1 and Z2 is the entire space X.
    union_of_subsets = Z1.union(Z2)

    print("We can find two proper closed subsets, Z1 and Z2:")
    print(f"Z1 = {Z1}")
    print(f"Z2 = {Z2}")
    
    print("\nThe union of these subsets is X itself, demonstrating the space is not irreducible.")
    
    # To satisfy the "output each number in the final equation" prompt,
    # we will print the elements of the sets in the final equation.
    # We sort the elements for a consistent, readable output.
    set_x_str = ", ".join(sorted(list(X)))
    set_z1_str = ", ".join(sorted(list(Z1)))
    set_z2_str = ", ".join(sorted(list(Z2)))

    print("\nFinal Equation:")
    print(f"{{{set_x_str}}} = {{{set_z1_str}}} U {{{set_z2_str}}}")

    print("\n--- Conclusion ---")
    n = 2
    print(f"Since a 2-point space can be constructed to be not irreducible,")
    print("and 0- and 1-point spaces are always irreducible,")
    print(f"the smallest non-negative integer n is {n}.")


solve()
<<<2>>>