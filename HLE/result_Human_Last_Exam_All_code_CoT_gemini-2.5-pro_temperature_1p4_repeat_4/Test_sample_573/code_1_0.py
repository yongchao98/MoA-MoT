import math

def calculate_simplices():
    """
    Calculates the number of n-simplices for n <= 5, N=200, and k=13.

    An n-simplex in the nerve of the overcategory N_.(Z_N)_{k/} corresponds to a
    non-decreasing sequence of n+1 integers j_0, ..., j_n such that
    -N <= j_0 <= ... <= j_n <= k.

    This is a problem of combinations with repetition. The number of integers
    to choose from is M = k - (-N) + 1 = N + k + 1. We are choosing n+1 integers.
    The formula is C(M + (n+1) - 1, n+1) = C(N + k + n + 1, n+1).
    """
    N = 200
    k = 13

    # M is the number of integers in the range [-N, k]
    M = N + k + 1

    print(f"For N = {N} and k = {k}, the number of n-simplices is calculated as C({M} + n, n + 1):")
    print("-" * 60)

    for n in range(6):  # n from 0 to 5
        # The number of items to choose from in the binomial coefficient
        num = M + n
        # The number of items to choose
        den = n + 1

        # Calculate the number of n-simplices using the combinations formula
        try:
            num_simplices = math.comb(num, den)
            print(f"For n = {n}, the number of {n}-simplices is C({num}, {den}) = {num_simplices}")
        except ValueError as e:
            print(f"Could not calculate for n = {n}: {e}")

calculate_simplices()