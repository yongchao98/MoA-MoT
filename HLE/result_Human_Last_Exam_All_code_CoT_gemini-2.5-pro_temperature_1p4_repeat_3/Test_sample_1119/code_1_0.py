import math

def count_sequences(N, K, M):
    """
    Calculates the number of possible sequences based on the given constraints.

    The number of sequences is given by the formula:
    Sum_{j=0}^{K-1} (-1)^j * C(K-1, j) * C(N - M*j, K)
    
    where C(n, k) is the binomial coefficient "n choose k".
    """
    
    print(f"Calculating for N={N}, K={K}, M={M}")

    # Check the given condition M(K-1) < N
    if not (M * (K - 1) < N):
        print(f"Warning: The condition M(K-1) < N (i.e., {M*K-M} < {N}) is not met.")
        print("The formula is derived assuming this condition holds.")

    # A safe combination function that returns 0 if n < k
    def nCr_safe(n, r):
        if r < 0 or r > n:
            return 0
        return math.comb(n, r)

    terms = []
    total_count = 0
    
    # Calculate each term in the inclusion-exclusion sum
    for j in range(K):
        # Calculate C(K-1, j) * C(N - M*j, K)
        term_value = nCr_safe(K - 1, j) * nCr_safe(N - M * j, K)
        terms.append(term_value)
        
        # Add or subtract the term from the total
        if j % 2 == 0:  # for j = 0, 2, 4, ...
            total_count += term_value
        else:  # for j = 1, 3, 5, ...
            total_count -= term_value
    
    # --- Output the results ---
    print("\nThe number of sequences is the result of the following equation:")
    
    # Build the equation string, e.g., "15504 - 24752 + 12012 - 1844 + 56 = 976"
    equation_str = ""
    # Find the last non-zero term to avoid printing trailing "+ 0 - 0"
    last_nonzero_idx = -1
    for i in range(len(terms) - 1, -1, -1):
        if terms[i] != 0:
            last_nonzero_idx = i
            break
            
    if last_nonzero_idx == -1:
        equation_str = "0 = 0"
    else:
        # Build the string from the list of term values
        equation_parts = [str(terms[0])]
        for i in range(1, last_nonzero_idx + 1):
            if i % 2 == 1:
                equation_parts.append("-")
            else:
                equation_parts.append("+")
            equation_parts.append(str(terms[i]))
        equation_str = " ".join(equation_parts) + f" = {total_count}"

    print(equation_str)
    
    print(f"\nFinal Answer: The total number of possible sequences is {total_count}.")
    return total_count

if __name__ == '__main__':
    # Example values for N, K, and M
    # N: the maximum value for any number in the sequence.
    # K: the length of the sequence.
    # M: the maximum increase between consecutive numbers.
    N = 20
    K = 5
    M = 3
    
    count_sequences(N, K, M)
