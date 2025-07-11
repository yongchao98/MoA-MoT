def count_non_compact_roots():
    """
    Calculates the number of non-compact positive roots for a real form of C_8.

    The Vogan diagram W--B--W--B--B--W--B==B means that the simple roots
    alpha_1, alpha_3, and alpha_6 are "white". A positive root is non-compact
    if the sum of its coefficients corresponding to these white simple roots is odd.
    """
    n = 8
    # The indices of the white simple roots are 1, 3, 6.
    S_W = {1, 3, 6}
    
    # --- Count non-compact short roots (type e_i - e_j) ---
    # Root expression: sum_{k=i}^{j-1} alpha_k
    # All coefficients are 1. We need |S_W intersect [i, j-1]| to be odd.
    short_non_compact = 0
    for j in range(2, n + 1):
        for i in range(1, j):
            # Interval of indices is [i, j-1]
            intersection_size = 0
            for p in S_W:
                if i <= p < j:
                    intersection_size += 1
            if intersection_size % 2 == 1:
                short_non_compact += 1

    # --- Count non-compact long roots ---
    
    # Type 2e_i
    # Expression: 2*sum_{k=i}^{n-1} alpha_k + alpha_n (for i < n) or alpha_n (for i=n)
    # The coefficients of alpha_k for k in {1,3,6} are either 0 or 2, which are even.
    # Therefore, the sum c_1+c_3+c_6 is always even.
    long_non_compact_2e = 0

    # Type e_i + e_j
    long_non_compact_eiej = 0
    # Case 1: 1 <= i < j < n
    # Expression: sum_{k=i}^{j-1} alpha_k + 2*sum_{k=j}^{n-1} alpha_k + alpha_n
    # The sum of white coeffs is |S_W intersect [i, j-1]| + 2*|S_W intersect [j, n-1]|.
    # The parity depends only on |S_W intersect [i, j-1]| being odd.
    for j in range(2, n): # j from 2 to 7
        for i in range(1, j):
            intersection_size = 0
            for p in S_W:
                if i <= p < j:
                    intersection_size += 1
            if intersection_size % 2 == 1:
                long_non_compact_eiej += 1

    # Case 2: 1 <= i < j = n
    # Expression: sum_{k=i}^{n-1} alpha_k + alpha_n
    # The sum of white coeffs is |S_W intersect [i, n-1]|.
    for i in range(1, n): # i from 1 to 7
        intersection_size = 0
        for p in S_W:
            if i <= p < n: # i.e., p is in [i, 7]
                intersection_size += 1
        if intersection_size % 2 == 1:
            long_non_compact_eiej += 1

    total_non_compact = short_non_compact + long_non_compact_2e + long_non_compact_eiej
    
    print(f"Number of non-compact short roots (type e_i - e_j): {short_non_compact}")
    print(f"Number of non-compact long roots (type 2e_i): {long_non_compact_2e}")
    print(f"Number of non-compact long roots (type e_i + e_j): {long_non_compact_eiej}")
    print(f"Total number of non-compact positive roots = {short_non_compact} + {long_non_compact_2e} + {long_non_compact_eiej} = {total_non_compact}")

count_non_compact_roots()