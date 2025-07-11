import math

def solve():
    """
    Calculates the number of n-simplices for the given mathematical construction.
    """
    N = 200
    k = 13
    n_values = [0, 1, 2, 3, 4, 5]

    # The objects of the over-category Z_N,k/ are the integers m
    # such that k <= m <= N.
    # The size of this set of objects is M.
    M = N - k + 1

    print(f"Given N = {N} and k = {k}.")
    print(f"The set of objects in the over-category is {{m | {k} <= m <= {N}}}.")
    print(f"The size of this set is M = {N} - {k} + 1 = {M}.")
    print("The number of n-simplices is the number of non-decreasing sequences of length n+1.")
    print(f"The formula for this is C(M + n, n + 1) = C({M} + n, n + 1).\n")

    results = []
    for n in n_values:
        # The number of n-simplices is the number of ways to choose n+1 elements
        # from a set of size M with replacement. This is a combination with repetition problem.
        # The formula is C(M + (n+1) - 1, n+1) = C(M + n, n+1).
        comb_n = M + n
        comb_k = n + 1
        
        # Calculate the number of simplices
        num_simplices = math.comb(comb_n, comb_k)
        results.append(num_simplices)

        # Print the result including the equation
        print(f"For n = {n}, the number of {n}-simplices is C({comb_n}, {comb_k}) = {num_simplices}")

if __name__ == "__main__":
    solve()