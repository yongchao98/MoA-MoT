def solve_ring_problem():
    """
    Solves the problem by enumerating the known mathematical results
    and calculating the final count.
    """
    # The problem asks for the size of a set of rings which have unique length factorizations (Half-Factorial Domains or HFDs).
    # The set is the union of two types of rings for square-free d > 0:
    # 1. The ring of integers O(Q(sqrt(-d))), which is the maximal order.
    # 2. The ring Z[sqrt(-d)] when it's not the maximal order (i.e., when d = 3 mod 4).

    # Step 1: Count HFDs among maximal orders O(Q(sqrt(-d))).
    # For these rings (Dedekind domains), being an HFD is equivalent to having a class number h(-d) of 1 or 2.

    # These are the square-free d > 0 for which h(-d) = 1 (9 rings).
    d_h1 = [1, 2, 3, 7, 11, 19, 43, 67, 163]
    num_h1 = len(d_h1)
    print(f"For maximal orders, the HFD property is equivalent to class number 1 or 2.")
    print(f"Number of rings with class number 1: {num_h1}")
    print(f"The values of d are: {d_h1}\n")

    # These are the square-free d > 0 for which h(-d) = 2 (21 rings).
    d_h2 = [5, 6, 10, 13, 15, 17, 21, 22, 33, 35, 37, 51, 58, 91, 115, 123, 187, 235, 267, 403, 427]
    num_h2 = len(d_h2)
    print(f"Number of rings with class number 2: {num_h2}")
    print(f"The values of d are: {d_h2}\n")

    total_maximal_hfd = num_h1 + num_h2
    print(f"Total maximal orders that are HFDs = {num_h1} + {num_h2} = {total_maximal_hfd}\n")

    # Step 2: Count HFDs among the non-maximal orders Z[sqrt(-d)] (where d is square-free and d = 3 mod 4).
    # From known theorems, only two such rings are HFDs.
    d_non_maximal_hfd = [3, 7]
    num_non_maximal_hfd = len(d_non_maximal_hfd)
    print(f"Number of specified non-maximal orders that are HFDs: {num_non_maximal_hfd}")
    print(f"The values of d are: {d_non_maximal_hfd}\n")
    
    # Step 3: Calculate the total size of the union.
    # The sets of rings are disjoint. The maximal order O(Q(sqrt(-d))) is never equal
    # to the non-maximal order Z[sqrt(-d)] (for d=3, 7).
    # Thus, we simply add the counts.
    total_rings = total_maximal_hfd + num_non_maximal_hfd
    
    print("The sets of maximal and non-maximal HFDs are disjoint.")
    print("The final equation is the sum of the sizes of these sets:")
    print(f"{num_h1} (h=1) + {num_h2} (h=2) + {num_non_maximal_hfd} (non-maximal) = {total_rings}")

solve_ring_problem()