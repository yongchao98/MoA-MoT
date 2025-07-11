def count_non_compact_roots():
    """
    Calculates the number of non-compact positive roots for the given real form of C8.
    """
    n = 8
    white_indices = {1, 3, 6}
    
    # --- Type 1 ---
    # Roots: sum_{k=i}^{j-1} alpha_k for 1 <= i < j <= 8
    count1 = 0
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            sum_white_coeffs = 0
            for k in range(i, j):
                if k in white_indices:
                    sum_white_coeffs += 1
            if sum_white_coeffs % 2 == 1:
                count1 += 1

    # --- Type 2 ---
    # Roots: sum_{k=i}^{j-1} alpha_k + 2*sum_{k=j}^{7} alpha_k + alpha_8 for 1 <= i < j <= 7
    count2 = 0
    for i in range(1, n):
        for j in range(i + 1, n):
            sum_white_coeffs = 0
            # Parity only depends on coeffs that are 1
            for k in range(i, j):
                if k in white_indices:
                    sum_white_coeffs += 1
            if sum_white_coeffs % 2 == 1:
                count2 += 1
    
    # --- Type 3 ---
    # Roots: sum_{k=i}^{7} alpha_k + alpha_8 for 1 <= i <= 7
    count3 = 0
    for i in range(1, n):
        sum_white_coeffs = 0
        for k in range(i, n):
            if k in white_indices:
                sum_white_coeffs += 1
        if sum_white_coeffs % 2 == 1:
            count3 += 1

    # --- Type 4 ---
    # Roots: 2*sum_{k=i}^{7} alpha_k + alpha_8 for 1 <= i <= 7
    # Sum of white coeffs is always even.
    count4 = 0

    # --- Type 5 ---
    # Root: alpha_8
    # Sum of white coeffs is 0.
    count5 = 0
    
    total_non_compact = count1 + count2 + count3 + count4 + count5
    
    print(f"Number of non-compact roots of Type 1: {count1} out of 28")
    print(f"Number of non-compact roots of Type 2: {count2} out of 21")
    print(f"Number of non-compact roots of Type 3: {count3} out of 7")
    print(f"Number of non-compact roots of Type 4: {count4} out of 7")
    print(f"Number of non-compact roots of Type 5: {count5} out of 1")
    print(f"Total number of non-compact positive roots is: {count1} + {count2} + {count3} + {count4} + {count5} = {total_non_compact}")

count_non_compact_roots()