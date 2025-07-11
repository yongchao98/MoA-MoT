def solve_irreducible_space_problem():
    """
    This function explains and solves the problem of finding the smallest n
    for which an n-point topological space can be not irreducible.
    """
    n_answer = 2

    print("To find the smallest non-negative integer n for which an n-point topological space is not irreducible (i.e., is reducible), we'll analyze small values of n.")
    print("-" * 50)

    print("Case n = 0:")
    print("The space is X = ∅ (the empty set). The only closed subset is ∅ itself. There are no proper closed subsets, so it cannot be written as a union of them. Thus, any 0-point space is irreducible.")
    print("-" * 50)

    print("Case n = 1:")
    print("The space is X = {p} (a single point). The only proper closed subset is the empty set, ∅. Any finite union of ∅ results in ∅, which is not X. Thus, any 1-point space is irreducible.")
    print("-" * 50)

    print(f"Case n = {n_answer}:")
    print(f"Let's consider a {n_answer}-point space, X = {{0, 1}}.")
    print("We can define the discrete topology on X, where every subset is open. Consequently, every subset is also closed.")
    print("\nThe closed sets in this space are: ∅, {0}, {1}, and {0, 1}.")
    print("The proper closed subsets are those not equal to X, which are: ∅, {0}, and {1}.")
    print("\nA space is reducible if it can be written as a finite union of proper closed subsets. Let's select two proper closed subsets:")
    print("Z1 = {0}")
    print("Z2 = {1}")
    print("\nTheir union covers the entire space X:")
    
    # Printing the final equation with the numbers
    z1_set_str = "{0}"
    z2_set_str = "{1}"
    x_set_str = "{0, 1}"
    print(f"    {x_set_str} = {z1_set_str} U {z2_set_str}")

    print("\nThis demonstrates that a 2-point space can be reducible.")
    print(f"Since spaces with n=0 and n=1 are always irreducible, the smallest such integer is {n_answer}.")


solve_irreducible_space_problem()

print("\n<<<2>>>")