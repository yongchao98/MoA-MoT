import math

def solve_sequence_count():
    """
    Calculates the number of possible sequences based on N, K, and M.

    A sequence of K positive integers (a_1, a_2, ..., a_K) must satisfy:
    1. 1 <= a_i <= N
    2. a_1 < a_2 < ... < a_K (strictly increasing)
    3. a_{i+1} - a_i <= M (increase between consecutive numbers is at most M)
    It is given that M*(K-1) < N.

    The problem is solved by transforming it into a combinatorial problem of
    counting integer solutions to an equation with constraints, which is then
    solved using the Principle of Inclusion-Exclusion.
    The resulting formula is: Sum_{j=0 to K-1} [(-1)^j * C(K-1, j) * C(N - j*M, K)]
    where C(n, k) is the binomial coefficient "n choose k".
    """
    # Parameters for the problem
    # You can change these values to solve for a different case.
    N = 20
    K = 4
    M = 5

    print(f"Problem: Find the number of sequences given N={N}, K={K}, M={M}.")
    print("-" * 30)

    # Check the pre-condition, although the formula works regardless
    if M * (K - 1) >= N:
        print("Warning: The condition M*(K-1) < N is not met.")
        print("The formula is still applicable, but the problem premise is violated.")

    total_sequences = 0
    equation_terms = []

    # Apply the Inclusion-Exclusion formula
    # Sum over j from 0 to K-1
    for j in range(K):
        # Calculate the binomial coefficient C(N - j*M, K)
        # math.comb(n, k) returns 0 if k > n, which correctly handles cases
        # where the term should be zero.
        n_val = N - j * M
        comb_n_k = math.comb(n_val, K)
        
        # Calculate the binomial coefficient C(K-1, j)
        comb_k_1_j = math.comb(K-1, j)
        
        # Calculate the j-th term in the summation
        term = ((-1)**j) * comb_k_1_j * comb_n_k
        
        equation_terms.append(term)
        total_sequences += term

    # Build the final equation string for printing
    if not equation_terms:
        equation_str = "0 = 0"
    else:
        # Start with the first term
        equation_str = str(equation_terms[0])
        # Add the subsequent terms with appropriate signs
        for term in equation_terms[1:]:
            if term >= 0:
                equation_str += f" + {term}"
            else:
                # Use " - " and the absolute value of the term
                equation_str += f" - {abs(term)}"
    
    # Add the final result
    equation_str += f" = {total_sequences}"

    print("The number of possible sequences is calculated as follows:")
    print(equation_str)


solve_sequence_count()
<<<1375>>>