import math

def solve_c8_noncompact_roots():
    """
    Calculates the number of non-compact positive roots for a real form of C8.
    The real form is specified by the Vogan diagram:
    W -- B -- W -- B -- B -- W -- B == B
    """
    n = 8
    
    # The Vogan diagram specifies which simple roots are "black".
    # 1-based indices of black nodes: {2, 4, 5, 7, 8}
    black_indices_1based = {2, 4, 5, 7, 8}
    # Convert to 0-based indices for list access
    black_indices_0based = {i - 1 for i in black_indices_1based}

    # A positive root is non-compact if the sum of its coefficients
    # for the black simple roots is odd.

    # --- Part 1: Roots of the form e_i - e_j ---
    # Formula: e_i - e_j = sum_{k=i}^{j-1} alpha_k, for 1 <= i < j <= n
    num_type1_roots = n * (n - 1) // 2
    non_compact_type1 = 0
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            # The coefficients are 1 for alpha_k where k is in [i, j-1].
            black_sum = 0
            for k_idx in range(i - 1, j - 1):
                if k_idx in black_indices_0based:
                    black_sum += 1
            
            if black_sum % 2 != 0:
                non_compact_type1 += 1

    # --- Part 2: Roots of the form e_i + e_j ---
    # Formula: e_i + e_j = (sum_{k=i}^{n-1} alpha_k) + (sum_{k=j}^{n-1} alpha_k) + alpha_n
    # for 1 <= i <= j <= n
    num_type2_roots = n * (n + 1) // 2
    non_compact_type2 = 0
    for i in range(1, n + 1):
        for j in range(i, n + 1):
            coeffs = [0] * n
            # Add coefficients from sum_{k=i}^{n-1} alpha_k
            for k in range(i, n):
                coeffs[k - 1] += 1
            # Add coefficients from sum_{k=j}^{n-1} alpha_k
            for k in range(j, n):
                coeffs[k - 1] += 1
            # Add coefficient for alpha_n
            coeffs[n - 1] += 1
            
            black_sum = 0
            for idx in black_indices_0based:
                black_sum += coeffs[idx]
            
            if black_sum % 2 != 0:
                non_compact_type2 += 1

    # --- Final Calculation and Output ---
    total_non_compact = non_compact_type1 + non_compact_type2
    
    print(f"The Lie algebra is of type C_8.")
    print(f"The black simple roots correspond to indices {sorted(list(black_indices_1based))}.")
    print("A positive root is non-compact if the sum of its coefficients for these black roots is odd.")
    print("")
    print(f"There are {num_type1_roots} positive roots of the form e_i - e_j.")
    print(f"Number of non-compact roots of this form: {non_compact_type1}")
    print("")
    print(f"There are {num_type2_roots} positive roots of the form e_i + e_j.")
    print(f"Number of non-compact roots of this form: {non_compact_type2}")
    print("")
    print(f"Total number of non-compact positive roots = {non_compact_type1} + {non_compact_type2} = {total_non_compact}")

solve_c8_noncompact_roots()