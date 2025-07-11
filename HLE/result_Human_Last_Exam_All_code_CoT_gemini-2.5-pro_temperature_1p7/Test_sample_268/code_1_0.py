import collections

def solve():
    """
    Calculates the number of non-compact positive roots for the real form of C_8
    given by the Vogan diagram W--B--W--B--B--W--B==B.
    """
    n = 8
    # The Vogan diagram specifies the black simple roots (1-based indices)
    black_indices = {2, 4, 5, 7, 8}

    # The positive roots of C_n in terms of simple roots are used:
    # Type 1: e_i - e_j for 1 <= i < j <= n
    # Type 2: e_i + e_j for 1 <= i < j <= n
    # Type 3: 2*e_i      for 1 <= i <= n

    # --- Type 1 roots: e_i - e_j ---
    # Coefficients: c_k = 1 for i <= k < j, and 0 otherwise.
    count_type1 = 0
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            coeffs = collections.defaultdict(int)
            for k in range(i, j):
                coeffs[k] = 1

            black_sum = 0
            for idx in black_indices:
                black_sum += coeffs[idx]

            if black_sum % 2 != 0:
                count_type1 += 1

    # --- Type 2 roots: e_i + e_j ---
    # Coefficients: c_k=1 for i<=k<j, c_k=2 for j<=k<n, c_n=1
    count_type2 = 0
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            coeffs = collections.defaultdict(int)
            for k in range(i, j):
                coeffs[k] += 1
            for k in range(j, n):
                coeffs[k] += 2
            coeffs[n] += 1

            black_sum = 0
            for idx in black_indices:
                black_sum += coeffs[idx]

            if black_sum % 2 != 0:
                count_type2 += 1

    # --- Type 3 roots: 2*e_i ---
    # Coefficients: c_k=2 for i<=k<n, c_n=1
    count_type3 = 0
    for i in range(1, n + 1):
        coeffs = collections.defaultdict(int)
        for k in range(i, n):
            coeffs[k] += 2
        coeffs[n] += 1

        black_sum = 0
        for idx in black_indices:
            black_sum += coeffs[idx]

        if black_sum % 2 != 0:
            count_type3 += 1

    total_non_compact = count_type1 + count_type2 + count_type3

    print(f"The number of non-compact positive roots is the sum of counts from each type:")
    print(f"Type 1 (e_i - e_j): {count_type1}")
    print(f"Type 2 (e_i + e_j): {count_type2}")
    print(f"Type 3 (2*e_i)    : {count_type3}")
    print(f"Total number of non-compact positive roots = {count_type1} + {count_type2} + {count_type3} = {total_non_compact}")

solve()