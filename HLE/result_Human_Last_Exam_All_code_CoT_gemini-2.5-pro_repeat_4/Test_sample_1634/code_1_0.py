def solve_reducible_space():
    """
    This function explains and demonstrates the solution to find the smallest n
    for which an n-point topological space is not irreducible.
    """
    print("Goal: Find the smallest non-negative integer n such that there exists an n-point topological space that is not irreducible.")
    print("-" * 80)
    print("A topological space X is NOT irreducible (i.e., it is reducible) if it can be written as a finite union of proper closed subsets.")
    print("Example: X = Z_1 U Z_2, where Z_1 and Z_2 are closed sets, and Z_1 != X, Z_2 != X.")

    print("\n--- Examining small values of n ---")

    print("\n[Case n=0]: X = {}. The only closed set is X itself. There are no proper closed subsets. Thus, the 0-point space is irreducible.")

    print("\n[Case n=1]: X = {p}. The only proper closed subset is the empty set, {}. The union of empty sets is not X. Thus, any 1-point space is irreducible.")

    print("\n[Case n=2]: Let's see if a 2-point space can be reducible.")
    print("Let the space be X = {a, b}.")
    print("We want to find proper closed subsets Z_1, Z_2 such that X = Z_1 U Z_2.")
    print("Let's try Z_1 = {a} and Z_2 = {b}. Their union is indeed {a, b} = X.")
    print("\nFor this to be a valid example, Z_1 and Z_2 must be closed sets.")
    print("We can achieve this with the discrete topology on X, where every subset is open.")
    print("In the discrete topology, the open sets are: {}, {a}, {b}, {a, b}.")
    print("A set is closed if its complement is open:")
    print(" - Complement of {a} is {b}, which is an open set. So, {a} is closed.")
    print(" - Complement of {b} is {a}, which is an open set. So, {b} is closed.")
    print("\nSo, in the 2-point space with the discrete topology, we found two proper closed subsets {a} and {b} whose union is the entire space X.")
    print("This space is therefore not irreducible.")

    print("\n--- Conclusion ---")
    print("The smallest n for which a non-irreducible space exists is 2.")
    print("The equation demonstrating this is:")

    # Define the sets for the final equation
    X_set = "{a, b}"
    Z1_set = "{a}"
    Z2_set = "{b}"

    # Print the final equation with all its components
    print(f"{X_set} = {Z1_set} U {Z2_set}")

solve_reducible_space()