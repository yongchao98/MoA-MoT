import math

def calculate_simplices():
    """
    Calculates the number of n-simplices for N=200, k=13, and n<=5.

    An n-simplex in N_•(Z_N)_k/ corresponds to a non-decreasing sequence of
    n+1 integers (x_0, ..., x_n) where k <= x_0 <= ... <= x_n <= N.

    This is a combinatorial problem of choosing n+1 numbers from the set
    {k, k+1, ..., N} with replacement.
    The size of this set is M = N - k + 1.
    The number of n-simplices is given by the multiset coefficient C(M + (n+1) - 1, n+1),
    which simplifies to C(M + n, n+1).
    """
    N = 200
    k = 13

    # The size of the set of choices {k, k+1, ..., N}
    M = N - k + 1

    print(f"Given N = {N} and k = {k}.")
    print(f"The number of available integers is M = {N} - {k} + 1 = {M}.")
    print("The number of n-simplices is C(M + n, n + 1).\n")

    results = []
    for n in range(6):
        num = M + n
        den = n + 1
        
        # Calculate the number of simplices using the combinations formula
        num_simplices = math.comb(num, den)
        results.append(num_simplices)

        print(f"For n = {n}:")
        print(f"The number of {n}-simplices is C({M} + {n}, {n} + 1) = C({num}, {den}) = {num_simplices}")

    return results

if __name__ == '__main__':
    calculate_simplices()
