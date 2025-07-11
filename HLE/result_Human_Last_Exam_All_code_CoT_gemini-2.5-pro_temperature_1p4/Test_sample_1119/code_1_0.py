import math

def count_sequences(N, K, M):
    """
    Calculates the number of possible sequences of K positive integers up to N.

    Args:
        N: The maximum value for any number in the sequence.
        K: The length of the sequence.
        M: The maximum increase between consecutive numbers.

    Returns:
        The total number of possible sequences.
    """
    if M * (K - 1) >= N:
        print("Warning: The condition M*(K-1) < N does not hold.")
        # The formula is still valid, but the problem context might change.

    total_sum = 0
    lhs_parts = []
    rhs_parts = []

    # The formula is Sum_{s=0}^{K-1} (-1)^s * C(K-1, s) * C(N - s*M, K)
    # We loop s from 0 up to K-1 (the loop range(K) covers this).
    for s in range(K):
        # The binomial coefficient C(n,k) is 0 if k > n.
        # math.comb handles this, returning 0 if the first argument is negative.
        comb1 = math.comb(K - 1, s)
        n_for_comb2 = N - s * M
        comb2 = math.comb(n_for_comb2, K)

        term_value = comb1 * comb2
        
        # We add non-zero terms to the equation string for clarity
        if term_value > 0:
            sign = ""
            if len(lhs_parts) > 0:  # Not the first term
                sign = " - " if s % 2 != 0 else " + "
            
            lhs_parts.append(f"{sign}C({K-1}, {s})*C({n_for_comb2}, {K})")
            rhs_parts.append(f"{sign}{comb1}*{comb2}")

        if s % 2 != 0:
            total_sum -= term_value
        else:
            total_sum += term_value

    equation_lhs = "".join(lhs_parts)
    equation_rhs = "".join(rhs_parts)
    
    print(f"For N={N}, K={K}, M={M}:")
    print("\nThe number of sequences is given by the equation:")
    print(f"{equation_lhs}")
    print(f"= {equation_rhs}")
    print(f"= {total_sum}")

    return total_sum

if __name__ == '__main__':
    # Example values for N, K, and M.
    # You can change these to your desired values.
    N = 20
    K = 5
    M = 3
    
    final_answer = count_sequences(N, K, M)
    # The final answer is also printed inside the function.
    # For automated extraction, you could print(f"<<<{final_answer}>>>") here.
