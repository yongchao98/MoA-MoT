def solve():
    """
    Calculates the number of true boolean expressions of exactly 5 symbols.
    The method partitions expressions by their lowest-precedence operator
    to avoid double counting.
    """

    # N3[L] stores a tuple (num_true, num_false) for expressions of length L
    # containing only T, F, !, () symbols.
    N3 = {}
    N3[1] = (1, 1)  # T, F
    N3[2] = (1, 1)  # !F, !T
    # For length 3: !!L -> (1,1), (L) -> (1,1). Total (2,2)
    N3[3] = (2, 2)

    # N23[L] stores counts for expressions of length L with no top-level '|'.
    N23 = {}
    N23[1] = N3[1]
    N23[2] = N3[2]
    # For length 3, add 'A&B' forms where A,B are from N3[1].
    # True: T&T (1). False: T&F, F&T, F&F (3).
    g2_len3_counts = (1, 3)
    N23[3] = (N3[3][0] + g2_len3_counts[0], N3[3][1] + g2_len3_counts[1]) # (2+1, 2+3) -> (3, 5)

    # 1. Count expressions of the form A|B (len(A)+len(B)=4)
    count_group1 = 0
    # Case: len(A)=1, len(B)=3. (and symmetric 3,1)
    ntA, nfA = N23[1]
    ntB, nfB = N23[3]
    total_A = ntA + nfA
    total_B = ntB + nfB
    # A|B is true unless both A and B are false.
    true_count_1_3 = total_A * total_B - nfA * nfB
    count_group1 += true_count_1_3 * 2  # For (1,3) and (3,1) splits

    # Case: len(A)=2, len(B)=2
    ntA, nfA = N23[2]
    ntB, nfB = N23[2]
    total_A = ntA + nfA
    total_B = ntB + nfB
    true_count_2_2 = total_A * total_B - nfA * nfB
    count_group1 += true_count_2_2

    # 2. Count expressions of the form A&B (len(A)+len(B)=4)
    count_group2 = 0
    # A,B must be from Group 3 type expressions.
    # Case: len(A)=1, len(B)=3. (and symmetric 3,1)
    ntA, nfA = N3[1]
    ntB, nfB = N3[3]
    # A&B is true only if both A and B are true.
    true_count_1_3 = ntA * ntB
    count_group2 += true_count_1_3 * 2 # For (1,3) and (3,1) splits

    # Case: len(A)=2, len(B)=2
    ntA, nfA = N3[2]
    ntB, nfB = N3[2]
    true_count_2_2 = ntA * ntB
    count_group2 += true_count_2_2

    # 3. Count expressions of length 5 using only !, (), T, F
    # The true expressions are: !!!!T, (!!T), !(!T), !!(T), ((T))
    count_group3 = 5

    total_count = count_group1 + count_group2 + count_group3
    
    print("The total number of true boolean expressions of length 5 is the sum of three disjoint groups:")
    print(f"1. Expressions with '|' as the root operator: {count_group1}")
    print(f"2. Expressions with '&' as the root operator (and no '|'): {count_group2}")
    print(f"3. Expressions with only '!', '()', and literals: {count_group3}")
    print(f"Total = {count_group1} + {count_group2} + {count_group3} = {total_count}")

solve()