def solve():
    """
    Calculates the number of non-compact positive roots for a real form of C_8
    given by the Vogan diagram W--B--W--B--B--W--B==B.
    """
    n = 8
    # Vogan Diagram: W-B-W-B-B-W-B==B
    # The indices of the non-compact simple roots (1-based).
    non_compact_indices = {2, 4, 5, 7, 8}

    # Create an indicator list `d` where d[k]=1 if alpha_k is non-compact, else 0.
    # We use 1-based indexing for convenience, so the list has size n+1.
    d = [0] * (n + 1)
    for i in non_compact_indices:
        d[i] = 1

    count_cat1 = 0
    count_cat2 = 0
    count_cat3 = 0

    # Category 1: roots of the form e_i - e_j for 1 <= i < j <= n
    # These correspond to the sum of simple roots from alpha_i to alpha_{j-1}.
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            coeffs = [0] * (n + 1)
            for k in range(i, j):
                coeffs[k] = 1
            
            non_compact_sum = sum(coeffs[k] * d[k] for k in range(1, n + 1))
            if non_compact_sum % 2 == 1:
                count_cat1 += 1

    # Category 2: roots of the form e_i + e_j for 1 <= i < j <= n
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            coeffs = [0] * (n + 1)
            if j < n:
                # Root is sum_{k=i to j-1} alpha_k + 2 * sum_{k=j to n-1} alpha_k + alpha_n
                for k in range(i, j):
                    coeffs[k] = 1
                for k in range(j, n):
                    coeffs[k] = 2
                coeffs[n] = 1
            else: # j == n
                # Root is sum_{k=i to n-1} alpha_k + alpha_n
                for k in range(i, n):
                    coeffs[k] = 1
                coeffs[n] = 1
            
            non_compact_sum = sum(coeffs[k] * d[k] for k in range(1, n + 1))
            if non_compact_sum % 2 == 1:
                count_cat2 += 1

    # Category 3: roots of the form 2*e_i for 1 <= i <= n
    for i in range(1, n + 1):
        coeffs = [0] * (n + 1)
        if i < n:
            # Root is 2 * sum_{k=i to n-1} alpha_k + alpha_n
            for k in range(i, n):
                coeffs[k] = 2
            coeffs[n] = 1
        else: # i == n
            # Root is alpha_n
            coeffs[n] = 1
        
        non_compact_sum = sum(coeffs[k] * d[k] for k in range(1, n + 1))
        if non_compact_sum % 2 == 1:
            count_cat3 += 1
    
    total_count = count_cat1 + count_cat2 + count_cat3
    print(f"Number of non-compact positive roots of type (e_i - e_j): {count_cat1}")
    print(f"Number of non-compact positive roots of type (e_i + e_j): {count_cat2}")
    print(f"Number of non-compact positive roots of type (2e_i): {count_cat3}")
    print(f"Total number of non-compact positive roots is {count_cat1} + {count_cat2} + {count_cat3} = {total_count}")
    print(f"\n<<<37>>>")

solve()