import math

def solve_simplices_count():
    """
    Calculates the number of n-simplices for the specified category.

    The problem asks for the number of n-simplices in the nerve of the
    overcategory N_dot(Z_N)_{k/}. An n-simplex in this context is a
    sequence of integers (x_0, x_1, ..., x_n) such that
    -N <= x_0 <= x_1 <= ... <= x_n <= k.

    This is a combination with repetition problem. We are choosing n+1 numbers
    from the set of integers {-N, ..., k}. The size of this set is M = k - (-N) + 1.
    The number of ways to do this is given by the multiset coefficient formula:
    C(M + (n+1) - 1, n+1) = C(M + n, n+1).
    """
    N = 200
    k = 13
    n_values = range(6)  # n <= 5, so n = 0, 1, 2, 3, 4, 5

    # The set of objects is the integers from -N to k.
    # The size of this set is M.
    M = k - (-N) + 1

    print(f"Given N = {N} and k = {k}.")
    print(f"The objects of the overcategory are integers from {-N} to {k}.")
    print(f"The number of available objects is M = {k} - ({-N}) + 1 = {M}.\n")

    for n in n_values:
        # The number of n-simplices is the number of ways to choose n+1 elements
        # from M options with replacement, which is C(M + n, n + 1).
        num_simplices = math.comb(M + n, n + 1)
        
        print(f"For n = {n}, the number of {n}-simplices is C({M} + {n}, {n} + 1) = C({M + n}, {n + 1}) = {num_simplices}")

solve_simplices_count()