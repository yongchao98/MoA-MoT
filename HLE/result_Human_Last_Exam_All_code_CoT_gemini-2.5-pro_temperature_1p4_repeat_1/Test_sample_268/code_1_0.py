def solve_c8_vogan():
    """
    Calculates the number of non-compact positive roots for a C8 Lie algebra
    with a given Vogan diagram.
    """
    n = 8
    # The Vogan diagram W--B--W--B--B--W--B==B defines which simple roots are non-compact.
    # Simple roots are alpha_1, ..., alpha_8.
    # W (compact): alpha_1, alpha_3, alpha_6
    # B (non-compact): alpha_2, alpha_4, alpha_5, alpha_7, alpha_8
    
    # We create a boolean list based on 1-based indexing for simple roots.
    # is_black[k-1] is True if alpha_k is black (non-compact).
    is_black = [
        False,  # alpha_1 (W)
        True,   # alpha_2 (B)
        False,  # alpha_3 (W)
        True,   # alpha_4 (B)
        True,   # alpha_5 (B)
        False,  # alpha_6 (W)
        True,   # alpha_7 (B)
        True    # alpha_8 (B)
    ]

    # Initialize counts for different types of positive roots.
    count1 = 0  # Type e_i - e_j
    count2a = 0 # Type e_i + e_j, j < n
    count2b = 0 # Type e_i + e_n
    count3 = 0  # Type 2*e_i

    # A root is non-compact if the sum of its coefficients corresponding to black simple roots is odd.

    # Category 1: Roots of the form e_i - e_j for 1 <= i < j <= 8
    # Simple root expansion: sum_{k=i}^{j-1} alpha_k. All coefficients are 1.
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            num_black_coeffs = sum(1 for k in range(i, j) if is_black[k - 1])
            if num_black_coeffs % 2 == 1:
                count1 += 1

    # Category 2a: Roots of the form e_i + e_j for 1 <= i < j < 8
    # Expansion: sum_{k=i}^{j-1} alpha_k + 2*sum_{k=j}^{n-1} alpha_k + alpha_n
    # The parity of the sum of black coefficients is determined by sum over [i,j-1] and alpha_n.
    # Since alpha_8 is black, we check if the number of black roots in [i, j-1] is even.
    for i in range(1, n):
        for j in range(i + 1, n):
            num_black_coeffs_part1 = sum(1 for k in range(i, j) if is_black[k - 1])
            if num_black_coeffs_part1 % 2 == 0:
                count2a += 1

    # Category 2b: Roots of the form e_i + e_n for 1 <= i < n
    # Expansion: sum_{k=i}^{n-1} alpha_k + alpha_n
    # Since alpha_8 is black, we check if the number of black roots in [i, n-1] is even.
    for i in range(1, n):
        num_black_coeffs = sum(1 for k in range(i, n) if is_black[k - 1])
        if num_black_coeffs % 2 == 0:
            count2b += 1

    # Category 3: Roots of the form 2*e_i for 1 <= i <= n
    # For i < n, expansion is 2*sum_{k=i}^{n-1} alpha_k + alpha_n. Sum of black coeffs is always odd.
    # For i = n, expansion is alpha_n. Sum of black coeffs is 1 (odd).
    # Therefore, all n roots of this type are non-compact.
    count3 = n

    total_non_compact = count1 + count2a + count2b + count3
    
    print("Decomposition of the non-compact positive root count:")
    print(f"Type e_i - e_j: {count1}")
    print(f"Type e_i + e_j (i < j < 8): {count2a}")
    print(f"Type e_i + e_8 (i < 8): {count2b}")
    print(f"Type 2*e_i: {count3}")
    print(f"\nTotal non-compact positive roots = {count1} + {count2a} + {count2b} + {count3} = {total_non_compact}")

solve_c8_vogan()