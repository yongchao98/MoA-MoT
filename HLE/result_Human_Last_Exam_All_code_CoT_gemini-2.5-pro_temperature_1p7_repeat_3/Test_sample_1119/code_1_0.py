def solve_sequence_count(N, K, M):
    """
    Calculates the number of sequences of K positive integers up to N with specified constraints.
    
    The sequence a_1, ..., a_K must satisfy:
    1. 1 <= a_1 < a_2 < ... < a_K <= N
    2. a_{i+1} - a_i <= M for i = 1, ..., K-1
    """

    def combinations(n, k):
        """
        Calculates the binomial coefficient C(n, k), also known as "n choose k".
        Returns 0 if k > n or k < 0.
        """
        if k < 0 or k > n:
            return 0
        if k == 0 or k == n:
            return 1
        # Choose the smaller of k and n-k for efficiency
        if k > n // 2:
            k = n - k
        
        # Iteratively compute C(n,k) = n * (n-1) * ... * (n-k+1) / k!
        # to avoid large intermediate numbers and floating point errors.
        res = 1
        for i in range(k):
            res = res * (n - i) // (i + 1)
        return res

    print(f"For N={N}, K={K}, M={M}:")
    # The condition M(K-1) < N is given. Let's verify it.
    if not (M * (K - 1) < N):
        print(f"Note: The given condition M*(K-1) < N, which is {M}*({K}-1) < {N}, is not met by the inputs.")
        print("-" * 20)

    total_sequences = 0
    calculation_parts = []
    
    # The formula is: Sum_{p=0}^{K-1} (-1)^p * C(K-1, p) * C(N - p*M, K)
    for p in range(K):
        sign = (-1)**p
        
        comb1_val = combinations(K - 1, p)
        
        # Argument for the second combination
        n_for_comb2 = N - p * M
        comb2_val = combinations(n_for_comb2, K)
        
        term_value = sign * comb1_val * comb2_val
        total_sequences += term_value

        # Skip adding terms that are zero to the final string representation
        if term_value == 0:
            continue

        comb1_str = f"C({K-1}, {p})"
        comb2_str = f"C({n_for_comb2}, {K})"
        
        if not calculation_parts: # This is the first term in the expression
            if sign == 1:
                calculation_parts.append(f"{comb1_str} * {comb2_str}")
            else: # sign is -1
                calculation_parts.append(f"-{comb1_str} * {comb2_str}")
        else: # These are subsequent terms
            if sign == 1:
                calculation_parts.append(f" + {comb1_str} * {comb2_str}")
            else: # sign is -1
                calculation_parts.append(f" - {comb1_str} * {comb2_str}")

    final_calculation_str = "".join(calculation_parts)
    # If all terms were zero, the result is 0
    if not final_calculation_str:
        final_calculation_str = "0"
        
    print("The number of possible sequences is calculated as:")
    print(f"{final_calculation_str} = {total_sequences}")


# --- Example Usage ---
# You can change these values to solve for different inputs.
# These values correspond to the worked example in the thought process.
N_val = 10
K_val = 4
M_val = 3

solve_sequence_count(N_val, K_val, M_val)
<<<108>>>