import math

def solve_sequence_count(N, K, M):
    """
    Calculates the number of possible integer sequences based on the problem's constraints.

    The number of sequences is found using a combinatorial formula derived from the
    Principle of Inclusion-Exclusion. The formula is:
    Sum_{t=0 to K-1} [(-1)^t * C(K-1, t) * C(N - t*M, K)]

    This function calculates the result and prints a detailed breakdown of the
    equation and its components.
    
    Args:
        N: The maximum value for any number in the sequence.
        K: The length of the sequence.
        M: The maximum increase between consecutive numbers.
    """

    print(f"Calculating the number of sequences for N={N}, K={K}, M={M}\n")

    def combinations(n, k):
        # C(n,k) is 0 if k > n. math.comb would raise an error.
        if k < 0 or k > n:
            return 0
        return math.comb(n, k)

    # The sum can be stopped when N - t*M < K.
    # t*M <= N - K  => t <= (N - K) / M
    # The upper limit for t in the summation.
    if M <= 0:
        # If M is not positive, and we need a_{i+1} - a_i >= 1, no solution exists.
        limit = -1
    else:
        limit = (N - K) // M

    total_sequences = 0
    numeric_eq_parts = []
    value_eq_parts = []
    
    # Iterate through t for the inclusion-exclusion formula
    for t in range(min(K - 1, limit) + 1):
        sign_char = "-" if t % 2 != 0 else "+"

        # Part 1: C(K-1, t)
        comb1 = combinations(K - 1, t)

        # Part 2: C(N - t*M, K)
        comb2_n = N - t * M
        comb2 = combinations(comb2_n, K)
        
        # Calculate the full term for this value of t
        term_value = comb1 * comb2
        
        if t > 0:
            numeric_eq_parts.append(f" {sign_char} C({K-1}, {t}) * C({comb2_n}, {K})")
            value_eq_parts.append(f" {sign_char} {comb1} * {comb2}")
        else:
            # First term does not have a leading sign
            numeric_eq_parts.append(f"C({K-1}, {t}) * C({comb2_n}, {K})")
            value_eq_parts.append(f"{comb1} * {comb2}")

        # Add or subtract from the total sum
        total_sequences += ((-1)**t) * term_value

    # Print the full equation breakdown
    if not numeric_eq_parts:
        print("No solutions possible with the given constraints.")
        print("\nFinal Answer: 0")
        return

    print("The number of sequences is calculated by the formula:")
    print("".join(numeric_eq_parts))
    
    print("\nWhich evaluates to:")
    print("= " + "".join(value_eq_parts))

    print("\nResulting in:")
    print(f"= {total_sequences}")

if __name__ == '__main__':
    # Example values for N, K, and M
    # You can change these values to solve for different inputs.
    N = 20
    K = 4
    M = 3

    # Validate the condition M*(K-1) < N
    if M * (K - 1) >= N:
        print(f"The condition M*(K-1) < N is not met for N={N}, K={K}, M={M}.")
        print("The formula may not be applicable.")
    else:
        solve_sequence_count(N, K, M)