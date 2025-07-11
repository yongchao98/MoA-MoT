import math

def solve_sequence_count():
    """
    Calculates the number of possible sequences based on N, K, and M.
    
    A sequence of K positive integers a_1, ..., a_K must satisfy:
    1. 1 <= a_i <= N
    2. a_1 < a_2 < ... < a_K
    3. a_{i+1} - a_i <= M
    4. M*(K-1) < N
    """
    
    # --- Parameters of the problem ---
    # N: The maximum value for any number in the sequence.
    # K: The length of the sequence.
    # M: The maximum allowed increase between consecutive numbers.
    # Please set the values for N, K, and M below.
    N = 20
    K = 5
    M = 3
    
    # --- Check the given condition ---
    if not (M * (K - 1) < N):
        print(f"The condition M*(K-1) < N is not met for N={N}, K={K}, M={M}.")
        print(f"Calculation: {M}*({K}-1) = {M*(K-1)}, which is not less than {N}.")
        return

    # --- Calculation using the inclusion-exclusion formula ---
    # Formula: Sum_{j=0 to K-1} (-1)^j * C(K-1, j) * C(N - j*M, K)
    
    total_sequences = 0
    
    # To display the calculation steps
    terms_str_list = []
    values_str_list = []
    
    # Loop for the sum in the formula
    for j in range(K):
        # Calculate binomial coefficients C(n, k) = n! / (k! * (n-k)!)
        # math.comb(n, k) returns 0 if k > n, which correctly handles terms
        # where the upper value in the combination becomes negative or too small.
        
        # C(K-1, j)
        comb1 = math.comb(K - 1, j)
        
        # C(N - j*M, K)
        n_val = N - j * M
        comb2 = math.comb(n_val, K)

        # Calculate the term in the sum
        term = ((-1)**j) * comb1 * comb2
        total_sequences += term
        
        # Prepare strings for formatted output
        sign = '+' if j % 2 == 0 else '-'
        terms_str_list.append(f"{sign} C({K-1}, {j})*C({n_val}, {K})")
        values_str_list.append(f"{sign} {comb1}*{comb2}")

    # --- Print the results ---
    print("The number of possible sequences is found using the formula derived from the Principle of Inclusion-Exclusion:")
    print("Result = Sum_{j=0 to K-1} [(-1)^j * C(K-1, j) * C(N - j*M, K)]")
    print(f"\nFor N={N}, K={K}, M={M}, the calculation is:")
    
    # Print the formula with values
    print(" ".join(terms_str_list).lstrip('+ '))
    
    # Print the formula with combination values
    print("= " + " ".join(values_str_list).lstrip('+ '))
    
    # Print the final result
    print(f"= {total_sequences}")
    
    # Return the final answer in the required format
    return total_sequences

if __name__ == '__main__':
    final_answer = solve_sequence_count()
    if final_answer is not None:
        print(f"\n<<< {final_answer} >>>")