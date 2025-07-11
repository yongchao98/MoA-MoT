def solve():
    """
    Calculates av_n^k(1324) for n=333 and k=3.

    Let av_n^k(p) denote the number of p-avoiding permutations of length n with k inversions.
    The problem is to determine av_333^3(1324).
    """

    n = 333
    k = 3
    pattern_1 = "1324"
    pattern_2 = "1234"

    # Step 1: Use the known equality av_n^k(1324) = av_n^k(1234).
    # This result from Claesson and Linusson (2011) means we can instead count
    # the number of 1234-avoiding permutations with k inversions.
    print(f"We need to find av_{n}^{k}({pattern_1}).")
    print(f"A known combinatorial result states that for any n and k:")
    print(f"av_n^k({pattern_1}) = av_n^k({pattern_2})")
    print(f"So, the problem is equivalent to finding av_{n}^{k}({pattern_2}).")
    print("-" * 20)

    # Step 2: Use the inequality LIS(pi) >= n - inv(pi).
    # LIS(pi) is the length of the Longest Increasing Subsequence of a permutation pi.
    # inv(pi) is the number of inversions of pi.
    print(f"For any permutation π of length n with k inversions, a theorem states:")
    print(f"LIS(π) >= n - k")
    print("-" * 20)

    # Step 3: Apply to our case n=333, k=3.
    min_lis_length = n - k
    print(f"For a permutation π of length n = {n} with k = {k} inversions:")
    print(f"LIS(π) >= {n} - {k} = {min_lis_length}")
    print("-" * 20)

    # Step 4: Check for the 1234 pattern.
    # A permutation contains the pattern 1234 if it has an increasing subsequence of length 4.
    print(f"Any permutation with an LIS of length {min_lis_length} contains an increasing subsequence of length 4, because {min_lis_length} >= 4.")
    print(f"This means that EVERY permutation of length {n} with {k} inversions must contain the {pattern_2} pattern.")
    print("-" * 20)

    # Step 5: Conclude the result.
    av_n_k_1234 = 0
    av_n_k_1324 = av_n_k_1234

    print(f"Therefore, the number of permutations of length {n} with {k} inversions that AVOID the {pattern_2} pattern is 0.")
    print(f"So, av_{n}^{k}({pattern_2}) = {av_n_k_1234}.")
    print(f"From the equality in the first step, we conclude:")
    print(f"av_{n}^{k}({pattern_1}) = {av_n_k_1324}")

solve()