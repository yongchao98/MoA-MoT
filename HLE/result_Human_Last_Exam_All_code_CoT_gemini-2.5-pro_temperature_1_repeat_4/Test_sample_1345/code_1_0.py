def calculate_max_complex_zeros(N):
    """
    Calculates the maximal possible number of complex zeros for the given N.

    A zero k is complex if Re(k_j) != 0 and Im(k_j) != 0 for all j.
    The formula for the maximal number of such zeros is:
    - 0 for N = 1
    - N * 2^(N-1) for N > 1
    """
    if not isinstance(N, int) or N <= 0:
        print(f"Error: N must be a positive integer. Got: {N}")
        return

    print(f"Calculating for N = {N}:")

    if N == 1:
        result = 0
        print("For N=1, the degree of the characteristic polynomial is odd, requiring at least one real root.")
        print("Therefore, not all solutions can be complex.")
        print(f"The maximal number of complex zeros is {result}.")
    else: # N > 1
        # The formula for the total number of solutions is N * 2^(N-1)
        term_N = N
        term_2 = 2
        exponent = N - 1
        power_of_2 = term_2 ** exponent
        result = term_N * power_of_2

        print("For N>1, the degree of the characteristic polynomial is even.")
        print("It is possible to choose a matrix A such that all solutions are complex.")
        print("The formula for the maximal number is N * 2**(N-1).")
        print(f"Calculation: {term_N} * {term_2}**({N}-1) = {term_N} * {power_of_2} = {result}")
        print(f"The maximal number of complex zeros is {result}.")

    print("-" * 30)


# --- Main execution ---
# Demonstrate the calculation for several values of N.
calculate_max_complex_zeros(1)
calculate_max_complex_zeros(2)
calculate_max_complex_zeros(3)
calculate_max_complex_zeros(4)
calculate_max_complex_zeros(5)