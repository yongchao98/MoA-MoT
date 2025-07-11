import sys

def solve():
    """
    Calculates the number of imaginary quadratic rings (maximal or specific non-maximal)
    that have unique lengths of factorization (i.e., are Half-Factorial Domains).
    """

    # Part 1: Count maximal orders (rings of integers) that are HFDs.
    # A ring of integers O_K for K = Q(sqrt(-d)) is an HFD if its class number is 1 or 2.
    # The lists of d for which this holds are known.

    # Square-free d > 0 for which the class number of Q(sqrt(-d)) is 1.
    d_h1 = [1, 2, 3, 7, 11, 19, 43, 67, 163]
    # Square-free d > 0 for which the class number of Q(sqrt(-d)) is 2.
    d_h2 = [5, 6, 10, 13, 15, 22, 35, 37, 51, 58, 91, 115, 123, 187, 235, 267, 403, 427]

    count_maximal_hfd = len(d_h1) + len(d_h2)

    # Part 2: Count non-maximal orders Z[sqrt(-d)] that are HFDs.
    # These rings are non-maximal for square-free d > 0 where d = 3 (mod 4).
    # The order is an HFD if its class number is <= 2.
    # The class number h(Z[sqrt(-d)]) relates to the field's class number h(Q(sqrt(-d))).
    # The order is an HFD only for the following d from the lists above:
    # 1. d = 3
    # 2. d > 3, d = 7 (mod 8)

    non_maximal_hfd_d_values = []

    # Case d=3 is a known HFD.
    if 3 in d_h1:
        non_maximal_hfd_d_values.append(3)

    # Find d > 3 and d = 7 (mod 8) from the known lists for h=1 and h=2.
    all_d_h_le_2 = d_h1 + d_h2
    for d in all_d_h_le_2:
        if d > 3 and d % 4 == 3 and d % 8 == 7:
            non_maximal_hfd_d_values.append(d)
    
    count_non_maximal_hfd = len(non_maximal_hfd_d_values)

    # The total is the sum of the counts from the two disjoint sets of rings.
    total_count = count_maximal_hfd + count_non_maximal_hfd

    # Print the final equation as requested.
    print(f"{count_maximal_hfd} + {count_non_maximal_hfd} = {total_count}")

solve()