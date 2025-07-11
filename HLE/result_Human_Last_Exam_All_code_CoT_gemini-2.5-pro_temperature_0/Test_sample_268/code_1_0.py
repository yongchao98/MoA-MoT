def solve():
    """
    Calculates the number of non-compact positive roots for a real form of C8
    given by a Vogan diagram.
    """
    n = 8
    # The Vogan diagram W--B--W--B--B--W--B==B corresponds to painted simple roots
    # alpha_2, alpha_4, alpha_5, alpha_7, alpha_8.
    # Indices are 1-based, so we use a set of integers.
    painted_indices = {2, 4, 5, 7, 8}

    num_non_compact_type1 = 0
    # Type 1 roots: sum_{k=i}^{j-1} alpha_k, for 1 <= i < j <= 8
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            coeffs = {k: 1 for k in range(i, j)}
            
            # Check non-compactness condition
            painted_coeff_sum = sum(coeffs.get(k, 0) for k in painted_indices)
            if painted_coeff_sum % 2 == 1:
                num_non_compact_type1 += 1

    num_non_compact_type2 = 0
    # Type 2 roots: sum_{k=i}^{j-1} alpha_k + 2*sum_{k=j}^{n-1} alpha_k + alpha_n, for 1 <= i < j <= 8
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            coeffs = {}
            # Part 1: sum_{k=i}^{j-1} alpha_k
            for k in range(i, j):
                coeffs[k] = 1
            # Part 2: 2*sum_{k=j}^{n-1} alpha_k
            for k in range(j, n):
                coeffs[k] = 2
            # Part 3: alpha_n
            coeffs[n] = 1
            
            # Check non-compactness condition
            painted_coeff_sum = sum(coeffs.get(k, 0) for k in painted_indices)
            if painted_coeff_sum % 2 == 1:
                num_non_compact_type2 += 1

    num_non_compact_type3 = 0
    # Type 3 roots: 2*sum_{k=i}^{n-1} alpha_k + alpha_n, for 1 <= i <= 8
    for i in range(1, n + 1):
        coeffs = {}
        # Part 1: 2*sum_{k=i}^{n-1} alpha_k
        for k in range(i, n):
            coeffs[k] = 2
        # Part 2: alpha_n
        coeffs[n] = 1

        # Check non-compactness condition
        painted_coeff_sum = sum(coeffs.get(k, 0) for k in painted_indices)
        if painted_coeff_sum % 2 == 1:
            num_non_compact_type3 += 1
            
    total_non_compact = num_non_compact_type1 + num_non_compact_type2 + num_non_compact_type3

    print(f"Number of non-compact roots of Type 1: {num_non_compact_type1}")
    print(f"Number of non-compact roots of Type 2: {num_non_compact_type2}")
    print(f"Number of non-compact roots of Type 3: {num_non_compact_type3}")
    print(f"Total number of non-compact positive roots = {num_non_compact_type1} + {num_non_compact_type2} + {num_non_compact_type3} = {total_non_compact}")

solve()