def solve_c8_noncompact_roots():
    """
    Calculates the number of non-compact positive roots for a real form of C8
    specified by a Vogan diagram.
    """
    n = 8
    # The Vogan diagram is W-B-W-B-B-W-B==B for alpha_1 to alpha_8.
    # We use 1-based indexing for simple roots, consistent with conventions.
    black_indices = {2, 4, 5, 7, 8}
    
    # We use 0-based indexing for lists in Python.
    black_indices_0based = {i - 1 for i in black_indices}

    num_noncompact_positive_roots = 0

    # There are n^2 = 64 positive roots in total for C_n.
    # They can be expressed in terms of simple roots alpha_k.

    # Type 1 roots: sum_{k=i}^{j-1} alpha_k, for 1 <= i < j <= n
    # These correspond to e_i - e_j.
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            coeffs = [0] * n
            for k in range(i, j):
                # k is 1-based here
                coeffs[k - 1] = 1
            
            sum_black_coeffs = 0
            for idx_0based in black_indices_0based:
                sum_black_coeffs += coeffs[idx_0based]
            
            if sum_black_coeffs % 2 != 0:
                num_noncompact_positive_roots += 1

    # Type 2 roots: sum_{k=i}^{j-1} alpha_k + 2*sum_{k=j}^{n-1} alpha_k + alpha_n
    # for 1 <= i < j <= n. These correspond to e_i + e_j.
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            coeffs = [0] * n
            for k in range(i, j):
                coeffs[k - 1] = 1
            for k in range(j, n):
                coeffs[k - 1] = 2
            coeffs[n - 1] = 1 # Coefficient for alpha_n is 1
            
            sum_black_coeffs = 0
            for idx_0based in black_indices_0based:
                sum_black_coeffs += coeffs[idx_0based]
            
            if sum_black_coeffs % 2 != 0:
                num_noncompact_positive_roots += 1

    # Type 3 roots: 2*sum_{k=i}^{n-1} alpha_k + alpha_n, for 1 <= i <= n
    # These correspond to 2*e_i.
    for i in range(1, n + 1):
        coeffs = [0] * n
        for k in range(i, n):
            coeffs[k - 1] = 2
        coeffs[n - 1] = 1 # Coefficient for alpha_n is 1
        
        sum_black_coeffs = 0
        for idx_0based in black_indices_0based:
            sum_black_coeffs += coeffs[idx_0based]
            
        if sum_black_coeffs % 2 != 0:
            num_noncompact_positive_roots += 1

    print(f"The number of non-compact positive roots is {num_noncompact_positive_roots}.")

solve_c8_noncompact_roots()