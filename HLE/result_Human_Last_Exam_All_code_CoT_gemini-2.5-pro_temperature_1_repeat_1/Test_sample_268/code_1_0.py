def solve():
    """
    Calculates the number of non-compact positive roots for a real form of C8.

    The real form is specified by a Vogan diagram, which determines the set
    of non-compact simple roots. A positive root is non-compact if the sum
    of its coefficients corresponding to non-compact simple roots is odd.
    """
    n = 8
    # The Vogan diagram W--B--W--B--B--W--B==B indicates that the simple roots
    # alpha_2, alpha_4, alpha_5, alpha_7, alpha_8 are non-compact (black).
    S_B = {2, 4, 5, 7, 8}

    # ----- Family 1: Roots of the form e_i - e_j for 1 <= i < j <= n -----
    # These roots are given by the expression: sum_{k=i}^{j-1} alpha_k.
    # All coefficients are 1.
    count1 = 0
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            s_B_sum = sum(1 for k in range(i, j) if k in S_B)
            if s_B_sum % 2 != 0:
                count1 += 1

    # ----- Family 2: Roots of the form e_i + e_j for 1 <= i < j <= n -----
    count2 = 0
    for i in range(1, n):  # i from 1 to n-1
        for j in range(i + 1, n + 1):  # j from i+1 to n
            s_B_sum = 0
            if j < n:
                # Root is sum_{k=i}^{j-1} alpha_k + 2*sum_{k=j}^{n-1} alpha_k + alpha_n
                # Coefficients are 1 for k in [i, j-1], 2 for k in [j, n-1], 1 for k=n.
                # We only need to consider odd coefficients for the parity of the sum.
                s_B_sum = sum(1 for k in range(i, j) if k in S_B)
                if n in S_B:
                    s_B_sum += 1
            else:  # j == n
                # Root is sum_{k=i}^{n-1} alpha_k + alpha_n
                # All coefficients are 1 for k in [i, n].
                s_B_sum = sum(1 for k in range(i, n + 1) if k in S_B)
            
            if s_B_sum % 2 != 0:
                count2 += 1

    # ----- Family 3: Roots of the form 2e_i for 1 <= i <= n -----
    count3 = 0
    for i in range(1, n + 1):
        # If i < n, root is 2*sum_{k=i}^{n-1} alpha_k + alpha_n.
        # If i = n, root is alpha_n.
        # In both cases, the only coefficient that can be odd is for alpha_n, which is 1.
        # The parity of the sum depends only on whether alpha_n is non-compact.
        s_B_sum = 0
        if n in S_B:
            s_B_sum += 1
        
        if s_B_sum % 2 != 0:
            count3 += 1

    total_count = count1 + count2 + count3
    print(f"Number of non-compact roots from family e_i - e_j: {count1}")
    print(f"Number of non-compact roots from family e_i + e_j (i<j): {count2}")
    print(f"Number of non-compact roots from family 2e_i: {count3}")
    print(f"Total number of non-compact positive roots is: {count1} + {count2} + {count3} = {total_count}")

solve()