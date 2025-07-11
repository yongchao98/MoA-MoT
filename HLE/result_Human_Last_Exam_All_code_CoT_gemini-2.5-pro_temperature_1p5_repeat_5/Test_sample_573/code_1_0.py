import math

def solve_simplices_count():
    """
    Calculates the number of n-simplices for the given problem.

    The number of n-simplices in the nerve of the overcategory Z_N_{k/} is
    given by the number of sequences (j_0, ..., j_n) such that
    k <= j_0 <= j_1 <= ... <= j_n <= N.

    This is a classic combinatorial problem whose solution is C(M + n, n+1),
    where M = N - k + 1 is the number of objects in the overcategory.
    The formula can be written as C(N - k + n + 1, n + 1).
    """
    # Parameters from the problem statement
    N = 200
    k = 13
    ns = [0, 1, 2, 3, 4, 5]

    print(f"Calculating the number of n-simplices for N={N}, k={k}, and n <= 5.")
    print("-----------------------------------------------------------------")
    
    # M is the number of objects in the overcategory (Z_N)_{k/}
    # This is also the number of 0-simplices.
    M = N - k + 1

    for n in ns:
        # Number of (n+1)-element multisets from a set of size M.
        # This is equivalent to C(M + (n+1) - 1, n+1) = C(M + n, n+1)
        p = n + 1
        num_simplices = math.comb(M + n, p)

        # Print the full equation for clarity
        print(f"For n = {n}:")
        print(f"Number of {n}-simplices = C({N} - {k} + {n} + 1, {n} + 1) = C({M + n}, {p}) = {num_simplices}")
        print()

solve_simplices_count()