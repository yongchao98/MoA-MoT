def solve():
    """
    This script solves the given topological puzzle by translating it into an
    algebraic problem and finding the lexicographically minimal solution.
    """
    print("Step 1: Understand the condition for a manifold to be 'full'.")
    print("For a closed, orientable 4-manifold, being 'full' is equivalent to its")
    print("Euler characteristic (chi), signature (sigma), and second Stiefel-Whitney class (w2) all being zero.")

    print("\nStep 2: Analyze the building blocks M(a,b) = M(a) x M(b).")
    print("For any a, b >= 0, M(a,b) has sigma=0 and w2=0.")
    print("Its Euler characteristic is chi = 4 * (1-a) * (1-b).")
    print("Therefore, M(a,b) is 'full' if and only if its chi=0, which requires a=1 or b=1.")

    print("\nStep 3: Set up the problem for the connected sum.")
    print("We need a collection of M(ai, bi) that are NOT full, so ai!=1 and bi!=1 for all i.")
    print("Their connected sum M must be full. The sigma and w2 of M are automatically zero.")
    print("The final condition is that the Euler characteristic of the sum must be zero:")
    print("SUM(chi_i) = SUM(4 * (1-ai)*(1-bi)) = 0.")
    print("This simplifies to the core equation: SUM((1-ai)*(1-bi)) = 0.")

    print("\nStep 4: Find the minimal solution.")
    print("The number of manifolds, l, must be minimal. l=1 is impossible since (1-a)(1-b) != 0 for a,b != 1.")
    print("So, the minimal l is 2. The equation becomes: (1-a1)(1-b1) + (1-a2)(1-b2) = 0.")

    print("\nStep 5: Search for the lexicographically smallest tuple (a1, b1, a2, b2).")
    print("To do this, we find the lexicographically smallest valid pair (a1, b1), which determines")
    print("the required value for the term from (a2, b2). We then find the smallest such (a2, b2).")

    # Smallest valid pair (a1, b1) with a,b >= 0, !=1 and a1 <= b1 is (0,0).
    a1, b1 = 0, 0
    c1 = (1 - a1) * (1 - b1)
    print(f"\nThe smallest valid pair (a1, b1) is ({a1}, {b1}). This gives the term c1 = {c1}.")

    c2_target = -c1
    print(f"This requires the second term c2 to be {c2_target}.")

    # Find smallest valid pair (a2, b2) with a2 <= b2 that gives c2_target = -1.
    # By inspection, this must be (0,2), since (1-0)*(1-2) = -1.
    a2, b2 = 0, 2
    c2 = (1 - a2) * (1 - b2)
    print(f"The smallest valid pair (a2, b2) producing c2={c2_target} is ({a2}, {b2}).")

    # The two pairs that form the solution are (0,0) and (0,2).
    # To form the lexicographically smallest tuple, we order the pairs.
    pair1 = (a1, b1)
    pair2 = (a2, b2)
    
    # Python's tuple comparison is lexicographical, so this check ensures correct order.
    if pair1 > pair2:
        pair1, pair2 = pair2, pair1

    final_tuple = pair1 + pair2

    # Unpack for final equation string
    t_a1, t_b1 = pair1
    t_a2, t_b2 = pair2
    t_c1 = (1 - t_a1) * (1 - t_b1)
    t_c2 = (1 - t_a2) * (1 - t_b2)

    print("\nFinal calculation verification:")
    print(f"(1 - {t_a1}) * (1 - {t_b1}) + (1 - {t_a2}) * (1 - {t_b2}) = {t_c1} + {t_c2} = 0")

    # The required format is a flat tuple with no spaces.
    flat_tuple_str = "(" + ",".join(map(str, final_tuple)) + ")"
    print(f"\nThe lexicographically least tuple is {flat_tuple_str}.")

solve()