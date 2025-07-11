def solve():
    """
    Calculates the number of 1324-avoiding permutations of length n=333 with k=3 inversions.
    For large n and small k, this number is independent of n.
    The method is based on constructing such permutations from smaller indecomposable 1324-avoiding blocks.
    """

    # c_k: number of "indecomposable" 1324-avoiding permutations with k inversions.
    # c_0: Identity permutation (length 1)
    # c_1: (2,1)
    c = {
        0: 1, # Base case, identity
        1: 1, # sigma=(2,1)
        2: 2, # sigma=(2,3,1), (3,1,2)
    }

    # For k=3:
    # From S_3: (3,2,1) is one such permutation.
    # From S_4: There are 6 perms with 3 inversions.
    # Check 1432 for 1324 pattern: 1<3<4<2, it contains the pattern.
    # The other 5 (2341, 2413, 3142, 3214, 4123) avoid it.
    # It can be shown that no indecomposable 1324-avoider with 3 inversions exists for length > 4.
    c[3] = 1 + 5

    k = 3
    result = 0

    # Case 1: One block of k inversions (at beginning or end)
    # sigma oplus id_...
    case1 = c[k]
    # id_... oplus sigma
    case2 = c[k]
    result += case1 + case2
    print(f"Number of permutations with all 3 inversions at the beginning: {case1}")
    print(f"Number of permutations with all 3 inversions at the end: {case2}")

    # Case 2: Two blocks of inversions (split k = i + (k-i))
    # sigma_1 oplus ... oplus sigma_2
    case3 = 0
    for i in range(1, k // 2 + 1):
        j = k - i
        if i == j:
            # If partitions are same size, say k=1+1, it's c_1 * c_1.
            # sigma_1 oplus id oplus sigma_2. The blocks are distinct in position.
            term = c[i] * c[j]
            case3 += term
        else:
            # Partitions i and j are different, e.g., k=1+2.
            # (c_i ... c_j) and (c_j ... c_i) are distinct compositions.
            term = 2 * c[i] * c[j]
            case3 += term
    
    # Manually breaking it down for k=3 (1+2 split)
    split_1_2 = c[1] * c[2]
    split_2_1 = c[2] * c[1]
    print(f"Number of permutations with a split of (1, 2) inversions: {split_1_2}")
    print(f"Number of permutations with a split of (2, 1) inversions: {split_2_1}")
    result += split_1_2 + split_2_1

    # Case 3: Three blocks of inversions (k = 1+1+1)
    split_1_1_1 = c[1] * c[1] * c[1]
    print(f"Number of permutations with a split of (1, 1, 1) inversions: {split_1_1_1}")
    result += split_1_1_1
    
    print("\nTotal number of permutations:")
    print(f"{case1} (beg) + {case2} (end) + {split_1_2} (1,2-split) + {split_2_1} (2,1-split) + {split_1_1_1} (1,1,1-split) = {result}")

solve()