import math

def solve_sequence_count():
    """
    Calculates the number of possible sequences based on N, K, and M.

    A sequence of K positive integers (a_1, a_2, ..., a_K) must satisfy:
    1. 1 <= a_1 < a_2 < ... < a_K <= N
    2. a_{i+1} - a_i <= M for i in 1..K-1
    3. The condition M*(K-1) < N holds.
    """
    # --- Parameters of the problem ---
    # You can change these values to solve for a different case.
    N = 10
    K = 3
    M = 4

    # --- Plan Explanation ---
    print("This script calculates the number of sequences of K positive integers up to N,")
    print("where each number is greater than the one before it, and the increase between")
    print("consecutive numbers does not exceed M.")
    print("\nThe number of sequences is given by the inclusion-exclusion formula:")
    print("Count = sum_{j=0 to floor((N-K)/M)} (-1)^j * C(K-1, j) * C(N - j*M, K)\n")
    print("-" * 40)
    print(f"Solving for N={N}, K={K}, M={M}")
    print("-" * 40)

    # --- Calculation ---
    total_sequences = 0
    # The sum runs as long as the argument of the second combination is non-negative.
    # N - j*M >= K  =>  N - K >= j*M  =>  j <= (N-K)/M
    limit = (N - K) // M

    # Store parts of the final equation string for a clean final output
    equation_parts = []

    for j in range(limit + 1):
        # Calculate the two binomial coefficients for the current term
        comb1 = math.comb(K - 1, j)
        
        n_for_comb2 = N - M * j
        # C(n, k) is 0 if n < k, so we can handle it this way
        if n_for_comb2 < K:
            comb2 = 0
        else:
            comb2 = math.comb(n_for_comb2, K)

        term = comb1 * comb2

        if j % 2 == 1:
            # Odd j corresponds to subtraction in the formula
            total_sequences -= term
            sign_char = "-"
            signed_term = -term
        else:
            # Even j corresponds to addition
            total_sequences += term
            sign_char = "+"
            signed_term = term

        # Output each number and operation in the formula
        print(f"j={j}: term = {sign_char} C({K-1}, {j}) * C({n_for_comb2}, {K})")
        print(f"          = {sign_char} {comb1} * {comb2} = {signed_term}")

        # Add the numeric term to our equation string
        if j == 0:
            equation_parts.append(f"{signed_term}")
        else:
            equation_parts.append(f"{sign_char} {abs(signed_term)}")

    print("-" * 40)
    # Join the parts to form a readable equation
    final_equation_str = " ".join(equation_parts)
    print(f"Total number of sequences = {final_equation_str} = {total_sequences}")


# Execute the function
solve_sequence_count()
<<<80>>>