import math

def solve_task(N):
    """
    Calculates the maximal possible number of complex zeros for the given matrix system.
    The formula is N * 2^(N-1).
    """
    if not isinstance(N, int) or N < 1:
        print("Error: N must be a positive integer.")
        return

    # For N=1, the number of complex zeros is 0. The formula gives 1 * 2^0 = 1.
    # The single root k1 = -A11 is real. So the number of complex zeros is 0.
    # The formula N * 2^(N-1) holds for N >= 2 for the total number of roots,
    # and represents the max number of complex roots. For N=1, the answer is 0.
    if N == 1:
        result = 0
        print("For N = 1, the equation is A11 + k1 = 0, which has one real root.")
        print("The maximal number of complex zeros is 0.")
        return

    # The formula for the maximal number of roots is N * 2^(N-1)
    term_N = N
    term_base = 2
    term_exponent = N - 1
    
    # Calculate the result
    result = term_N * (term_base ** term_exponent)

    # Output the explanation and the result
    print(f"For a system of size N = {N}:")
    print("The formula for the maximal number of complex zeros is N * 2^(N-1).")
    print(f"The components of the calculation are:")
    print(f"N = {term_N}")
    print(f"Base = {term_base}")
    print(f"Exponent = N - 1 = {term_exponent}")
    print(f"The maximal number of complex zeros is {term_N} * {term_base}^{term_exponent} = {result}")

# Example for N=4
N = 4
solve_task(N)
