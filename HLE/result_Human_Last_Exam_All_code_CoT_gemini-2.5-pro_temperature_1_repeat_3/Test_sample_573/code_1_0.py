import math

def calculate_simplices():
    """
    Calculates the number of n-simplices in the specified category nerve
    for n from 0 to 5, with N=200 and k=13.

    An n-simplex in N_.(Z_N)_k/ corresponds to a sequence of integers
    (x_0, x_1, ..., x_n) such that k <= x_0 <= x_1 <= ... <= x_n <= N.
    The number of such sequences is given by the combinatorial formula for
    multisets, which is C(M + K - 1, K), where M is the number of integers
    to choose from (N - k + 1) and K is the length of the sequence (n + 1).
    This simplifies to C(N - k + n + 1, n + 1).
    """
    N = 200
    k = 13

    print(f"Calculating the number of n-simplices for N={N} and k={k}, for n from 0 to 5.")
    print("-" * 70)

    for n in range(6):
        # The arguments for the binomial coefficient C(m, r)
        m = N - k + n + 1
        r = n + 1
        
        # Calculate the result using math.comb for precision with large numbers
        result = math.comb(m, r)
        
        # Print the detailed equation for each n
        print(f"For n={n}, the number of simplices is C({N} - {k} + {n} + 1, {n} + 1) = C({m}, {r}) = {result}")

calculate_simplices()