def count_non_compact_roots():
    """
    Calculates the number of non-compact positive roots for a real form of C8
    given by the Vogan diagram W--B--W--B--B--W--B==B.
    """
    n = 8
    # The Vogan diagram indicates which simple roots are non-compact (white).
    # Using 1-based indexing for alpha_1, ..., alpha_8.
    white_nodes = {1, 3, 6}

    total_non_compact_count = 0
    
    # Category 1: Roots of the form e_i - e_j for 1 <= i < j <= 8
    # There are n*(n-1)/2 = 28 such roots.
    # The expression is sum_{k=i}^{j-1} alpha_k.
    # Coefficients are 1 for k in [i, j-1] and 0 otherwise.
    num_type1_nc = 0
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            white_coeff_sum = 0
            for k in range(i, j): # Corresponds to interval [i, j-1]
                if k in white_nodes:
                    white_coeff_sum += 1
            if white_coeff_sum % 2 != 0:
                num_type1_nc += 1

    # Category 2: Roots of the form e_i + e_j for 1 <= i < j <= 8
    # There are n*(n-1)/2 = 28 such roots.
    # Expression: sum_{k=i}^{j-1} alpha_k + 2*sum_{k=j}^{n-1} alpha_k + alpha_n
    num_type2_nc = 0
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            white_coeff_sum = 0
            # Part 1: sum from k=i to j-1 has coefficient 1.
            for k in range(i, j):
                if k in white_nodes:
                    white_coeff_sum += 1
            # Part 2: sum from k=j to n-1 has coefficient 2.
            # A coefficient of 2 adds an even number to the sum, so it does not
            # change the parity. We can ignore these terms in the parity check.
            if white_coeff_sum % 2 != 0:
                num_type2_nc += 1
                
    # Category 3: Roots of the form 2*e_i for 1 <= i <= 8
    # There are n=8 such roots.
    # Expression: 2*sum_{k=i}^{n-1} alpha_k + alpha_n
    num_type3_nc = 0
    for i in range(1, n + 1):
        white_coeff_sum = 0
        # For i < n, the coefficients of alpha_k are 2.
        # This will always result in an even sum of white coefficients.
        # For i = n, the root is alpha_n. The sum of white coefficients is 0.
        if white_coeff_sum % 2 != 0:
            num_type3_nc += 1

    total_non_compact_count = num_type1_nc + num_type2_nc + num_type3_nc
    
    print("Counting non-compact positive roots:")
    print(f"Number from roots of type e_i - e_j: {num_type1_nc}")
    print(f"Number from roots of type e_i + e_j: {num_type2_nc}")
    print(f"Number from roots of type 2*e_i: {num_type3_nc}")
    print("\nFinal equation for the total number of non-compact positive roots:")
    print(f"{num_type1_nc} + {num_type2_nc} + {num_type3_nc} = {total_non_compact_count}")


count_non_compact_roots()
