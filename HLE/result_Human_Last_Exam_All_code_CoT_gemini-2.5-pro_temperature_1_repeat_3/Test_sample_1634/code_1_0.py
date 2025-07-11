def solve_irreducible_space_problem():
    """
    This script determines the smallest non-negative integer n such that
    an n-point topological space can be not irreducible (i.e., reducible).
    It does so by explaining the logic and demonstrating the case for n=2.
    """
    print("--- Problem ---")
    print("Find the smallest non-negative integer n such that there exists an n-point topological space that is not irreducible.")
    
    print("\n--- Definition ---")
    print("A topological space X is REDUCIBLE (not irreducible) if it can be written as a union of a finite number of proper closed subsets.")
    print("A subset Z is a 'proper closed subset' if Z is a closed set and Z is not equal to the whole space X.")
    print("This means X is reducible if X = Z_1 U Z_2 U ... U Z_k, where each Z_i is a closed set and Z_i != X.")

    print("\n--- Analysis of small n ---")

    print("\n[Case n = 0]")
    print("Let X = {}. The only closed set is {} itself. There are no proper closed subsets. Thus, the 0-point space is irreducible.")

    print("\n[Case n = 1]")
    print("Let X = {p}. The only proper closed subset is the empty set {}. Any union of empty sets is still {}, which does not equal X. Thus, any 1-point space is irreducible.")

    print("\n[Case n = 2]")
    print("Let X = {1, 2}. We will try to construct a topology that makes X reducible.")
    print("To make X reducible, we need to write it as a union of proper closed subsets. Let's try to find two such sets, Z_1 and Z_2.")
    print("The only way for the union of two proper subsets to be X = {1, 2} is if Z_1 = {1} and Z_2 = {2}.")
    print("This requires a topology where both {1} and {2} are closed sets.")
    print("If {1} is closed, its complement {2} must be open.")
    print("If {2} is closed, its complement {1} must be open.")
    print("This is satisfied by the discrete topology, where every subset is open (and therefore also closed).")

    print("\n--- Verification for n = 2 ---")
    X = {1, 2}
    Z1 = {1}
    Z2 = {2}

    print(f"Consider the space X = {X} with the discrete topology.")
    print(f"Let Z1 = {Z1} and Z2 = {Z2}.")
    print("In this topology, both Z1 and Z2 are closed sets.")
    print("Z1 is a proper subset of X, and Z2 is a proper subset of X.")
    
    union_of_subsets = set(Z1).union(Z2)

    # Format sets into strings for the equation output, as requested.
    z1_str = "{" + ", ".join(map(str, sorted(list(Z1)))) + "}"
    z2_str = "{" + ", ".join(map(str, sorted(list(Z2)))) + "}"
    union_str = "{" + ", ".join(map(str, sorted(list(union_of_subsets)))) + "}"

    print("\nWe check if their union equals X. The equation is:")
    print(f"{z1_str} U {z2_str} = {union_str}")

    if union_of_subsets == X:
        print("\nThe union of these two proper closed subsets is the entire space X.")
        print("Therefore, this 2-point space is reducible.")
    
    print("\n--- Conclusion ---")
    print("Since spaces with n=0 and n=1 are irreducible, and we have constructed a reducible space for n=2, the smallest such integer is 2.")

# Execute the function to print the explanation.
solve_irreducible_space_problem()