import math

def combinations(n, k):
    """
    Calculates the binomial coefficient "n choose k", C(n, k).
    Handles cases where n < k by returning 0.
    """
    if k < 0 or k > n:
        return 0
    # math.comb is efficient and available in Python 3.8+
    try:
        return math.comb(n, k)
    # Manual implementation for older Python versions
    except AttributeError:
        if k == 0:
            return 1
        if k > n // 2:
            k = n - k
        
        numerator = 1
        for i in range(k):
            numerator *= (n - i)
        
        denominator = math.factorial(k)
        return numerator // denominator

def solve_sequence_count(N, K, M):
    """
    Calculates the number of sequences based on the derived combinatorial formula.
    
    The number of sequences is given by:
    Sum_{j=0}^{K-1} (-1)^j * C(K-1, j) * C(N - j*M, K)
    
    This function prints the step-by-step evaluation of this formula.
    """
    print(f"Solving for N={N}, K={K}, M={M}")
    print(f"Condition M*(K-1) < N: {M}*({K}-1) = {M*(K-1)} < {N}, which is {M*(K-1) < N}.")
    print("-" * 30)
    
    total_sequences = 0
    equation_str_parts = []
    values_str_parts = []

    # Loop from j = 0 to K-1 to calculate the summation
    for j in range(K):
        # Calculate the two combination terms
        comb1 = combinations(K - 1, j)
        n_arg_for_comb2 = N - j * M
        comb2 = combinations(n_arg_for_comb2, K)
        
        # Calculate the full term value
        term_value = comb1 * comb2
        
        # Build the string representing the formula part, e.g., "+ C(3, 0) * C(20, 4)"
        sign_char = "-" if j % 2 != 0 else "+"
        if j == 0:
            # First term does not need a leading sign
            equation_str_parts.append(f"C({K-1}, {j}) * C({n_arg_for_comb2}, {K})")
        else:
            equation_str_parts.append(f" {sign_char} C({K-1}, {j}) * C({n_arg_for_comb2}, {K})")
            
        # Add the calculated value to the total sum with the correct sign
        if j % 2 == 0:
            total_sequences += term_value
        else:
            total_sequences -= term_value
        
        # Build the string representing the values, e.g., "+ 4845"
        if j == 0:
             values_str_parts.append(f"{term_value}")
        else:
             values_str_parts.append(f" {sign_char} {term_value}")

    # Print the detailed breakdown of the calculation
    print("The number of sequences is calculated using the formula:")
    print("".join(equation_str_parts))
    print(f"= {''.join(values_str_parts)}")
    print(f"= {total_sequences}")
    
# --- Main execution ---
if __name__ == "__main__":
    # You can change these values to solve for different parameters
    N = 20
    K = 4
    M = 5
    
    solve_sequence_count(N, K, M)
