def solve_vc_dimension():
    """
    Calculates and explains the VC dimension for the given first-order logic fragment.
    """
    # Number of unary predicates in the schema S
    num_predicates = 4

    print("--- Step 1: Understanding the Concept Class ---")
    print(f"The schema S has {num_predicates} unary predicates: P_1, P_2, P_3, P_4.")
    print("The logical fragment is FO_{exists, conjunction, top, bottom}[S].")
    print("A concept in this context is a formula phi(x) with one free variable x.")
    print("Any such formula phi(x) can be shown to be equivalent to one of two forms:")
    print("  1. The always-false concept (equivalent to bot).")
    print("  2. A monotone conjunction of the base predicates, like P_1(x) AND P_3(x).")
    print("So, the problem is to find the VC dimension of the class of monotone conjunctions over 4 boolean features.\n")

    # Let k be the number of predicates
    k = num_predicates

    print(f"--- Step 2: Proving VC dimension >= {k} ---")
    print(f"We will show that a set of {k} points can be shattered.")
    print("Let the points be {x_1, x_2, x_3, x_4}.")
    print("We can construct a model where for each predicate P_j, it is true for all points EXCEPT x_j.")
    print("Let's visualize this with a table (1 means True, 0 means False):")
    header = "       | " + " | ".join([f"P_{j+1}" for j in range(k)])
    print(header)
    print("-" * len(header))
    for i in range(k):
        # row[j] is 1 if P_{j+1}(x_{i+1}) is True
        row = [1 if i != j else 0 for j in range(k)]
        print(f"  x_{i+1}  |   " + "   |   ".join(map(str, row)))

    print("\nWith this setup, any subset Y of points can be defined by a formula.")
    print("The formula is the conjunction of all predicates P_j for which the corresponding point x_j is NOT in Y.")
    print("Example: To select only {x_1, x_4}, the points not in the set are {x_2, x_3}. The formula is P_2(x) AND P_3(x).")
    print(f"This construction can generate all 2^{k} = {2**k} possible labelings.")
    print(f"Therefore, a set of size {k} can be shattered, and the VC dimension is at least {k}.\n")

    # Let m be the size of the set we test for the upper bound
    m = k + 1
    print(f"--- Step 3: Proving VC dimension < {m} ---")
    print(f"We will show that no set of {m} points can be shattered, using a proof by contradiction.")
    print(f"Assume we can shatter a set of {m} points {x_1, ..., x_{m}}.")
    print("If the set is shattered, then for any point x_j, we must be able to label it as FALSE and all other points as TRUE.")
    print("Let's call the formula for this labeling phi_j(x).")
    print("Since phi_j(x) is a conjunction, for it to be FALSE at x_j, at least one predicate in it, say P_{l_j}, must be FALSE at x_j.")
    print("Since phi_j(x) must be TRUE for all other points x_i (where i != j), the predicate P_{l_j} must be TRUE for all those points.")
    print("\nSo, for each point x_j (from j=1 to 5), we've found a predicate index l_j (from 1 to 4) such that:")
    print(f"  - P_{{l_j}}(x_j) is FALSE")
    print(f"  - P_{{l_j}}(x_i) is TRUE for all i != j")
    
    print(f"\nThis gives us a mapping from {m} points to {k} predicate indices.")
    print(f"By the pigeonhole principle, since there are {m} points ('pigeons') and only {k} indices ('pigeonholes'), at least two points must map to the same index.")
    print("Let's say points x_a and x_b both map to the same predicate index, l*.")
    
    print("\nThis leads to a direct contradiction:")
    print("  - From the property for point x_a, P_{l*}(x_b) must be TRUE (because b is not a).")
    print("  - From the property for point x_b, P_{l*}(x_b) must be FALSE.")
    print("A predicate cannot be both TRUE and FALSE for the same point P_{l*}(x_b).")
    print(f"The contradiction proves our assumption was wrong. No set of {m} points can be shattered.")
    print(f"Therefore, the VC dimension must be less than {m}.\n")
    
    print("--- Step 4: Conclusion ---")
    vc_dimension = k
    print(f"From Step 2, we have: VC dimension >= {k}")
    print(f"From Step 3, we have: VC dimension < {k+1}")
    print(f"The only integer value that satisfies these conditions is {k}.")
    print("\nThe VC dimension of FO_{exists, and, top, bot}[S] with 4 unary predicates is:")
    # Final equation as requested by the user
    print(f"VC dimension = {vc_dimension}")

if __name__ == "__main__":
    solve_vc_dimension()