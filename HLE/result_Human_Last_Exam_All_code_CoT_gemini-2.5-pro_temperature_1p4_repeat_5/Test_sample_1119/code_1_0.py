import math

def count_sequences(N, K, M):
    """
    Calculates the number of possible sequences based on the given constraints.

    A sequence of K positive integers (a_1, ..., a_K) must satisfy:
    1. a_i <= N for all i.
    2. a_1 < a_2 < ... < a_K.
    3. a_{i+1} - a_i <= M for all i.
    4. The condition M*(K-1) < N holds.

    The number of such sequences is given by the formula:
    Sum_{j=0}^{K-1} [(-1)^j * C(K-1, j) * C(N - j*M, K)]
    where C(n, k) is "n choose k".
    """
    
    # Validate the condition given in the problem statement
    if not (M * (K - 1) < N):
        print(f"Warning: The condition M*(K-1) < N is not met. "
              f"M*(K-1) = {M*(K-1)}, N = {N}. The formula might still be applicable.")

    total_sequences = 0
    equation_parts = []

    # Iterate through j from 0 to K-1 as per the formula
    for j in range(K):
        # Calculate the components of the formula's term
        
        # First combination: C(K-1, j)
        comb_K_1_j = math.comb(K - 1, j)
        
        # Second combination: C(N - j*M, K)
        n_val = N - j * M
        
        # C(n, k) is 0 if n < k. math.comb handles this since Python 3.8.
        # Adding an explicit check for compatibility and clarity.
        if n_val < K:
            comb_njM_K = 0
        else:
            comb_njM_K = math.comb(n_val, K)
            
        # The value of the term for the sum (without the alternating sign yet)
        term_value = comb_K_1_j * comb_njM_K
        
        # Add the signed term to the total
        if j % 2 == 1: # Odd j corresponds to a negative term
            total_sequences -= term_value
        else: # Even j corresponds to a positive term
            total_sequences += term_value
        
        # Append the numeric value to our list for the final equation string
        equation_parts.append(str(term_value))

    # Build the final output string showing the calculation
    equation_str = " - ".join(equation_parts[:2])
    for i in range(2, len(equation_parts)):
        sign = "+" if i % 2 == 0 else "-"
        equation_str += f" {sign} {equation_parts[i]}"
        
    print("The number of possible sequences is given by the calculation:")
    print(f"{equation_str} = {total_sequences}")
    
    return total_sequences

if __name__ == '__main__':
    # Example values. You can change these to solve for a different case.
    N = 20
    K = 4
    M = 5
    
    print(f"Calculating for N={N}, K={K}, M={M}:")
    count_sequences(N, K, M)