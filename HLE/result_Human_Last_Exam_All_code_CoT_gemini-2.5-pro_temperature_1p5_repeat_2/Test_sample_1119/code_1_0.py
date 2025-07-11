import math

def solve_sequence_count(N, K, M):
    """
    Calculates the number of possible sequences based on the given constraints.

    A sequence of K positive integers (a_1, ..., a_K) must satisfy:
    1. 1 <= a_1 < a_2 < ... < a_K <= N
    2. a_{i+1} - a_i <= M for all i
    3. The condition M * (K - 1) < N holds.

    The number of such sequences is calculated using the formula:
    Sum_{j=0}^{K-1} [(-1)^j * C(K-1, j) * C(N - j*M, K)]
    """
    
    # Check the given condition for validity, although the formula works regardless.
    if not (M * (K - 1) < N):
        print(f"Warning: The condition M*(K-1) < N (i.e., {M}*({K}-1) < {N}) is not met.")
        print("The formula might result in 0 or an unexpected value for this case.")
        # return
        
    print(f"Calculating the number of sequences for N={N}, K={K}, M={M}\n")

    total_sequences = 0
    equation_str_parts = []
    value_str_parts = []

    # The summation is from j = 0 to K-1
    # We can technically loop until N - j*M < K, but math.comb handles it.
    for j in range(K): 
        # Calculate the binomial coefficients
        # C(n, k) is 0 if k > n, which handles cases where N - j*M < K
        try:
            comb1 = math.comb(K - 1, j)
            comb2 = math.comb(N - j * M, K)
        except ValueError:
            # math.comb throws ValueError for negative n, which we treat as 0 combinations.
            comb1 = math.comb(K - 1, j) if K-1 >= j else 0
            comb2 = 0

        # Determine the sign of the term
        sign = (-1) ** j
        term = sign * comb1 * comb2
        total_sequences += term

        # Build the string representation of the equation
        if j == 0:
            sign_char = ""
        elif sign > 0:
            sign_char = "+ "
        else:
            sign_char = "- "
            
        # Format strings for C(n,k)
        c1_str = f"C({K-1}, {j})"
        c2_str = f"C({N - j*M}, {K})"
        equation_str_parts.append(f"{sign_char}{c1_str} * {c2_str}")
        
        # Format string for values
        if j == 0:
             val_sign_char = "= "
        elif sign > 0:
             val_sign_char = "+ "
        else:
             val_sign_char = "- "
        
        value_str_parts.append(f"{val_sign_char}{comb1} * {comb2}")


    print("Formula:")
    print(' '.join(equation_str_parts))
    print("\nCalculation:")
    print(' '.join(value_str_parts))
    print(f"\nTotal number of possible sequences: {total_sequences}")


if __name__ == '__main__':
    # --- Example ---
    # You can change these values to solve for a different case.
    N = 20  # The numbers in the sequence cannot exceed N
    K = 5   # The sequence has K numbers
    M = 3   # The increase between consecutive numbers cannot exceed M

    # The condition M*(K-1) < N is 3*(5-1) = 12 < 20, which holds.
    
    solve_sequence_count(N, K, M)
    
    # --- Another Example (from the thought process) ---
    print("\n" + "="*40 + "\n")
    solve_sequence_count(N=10, K=3, M=4)