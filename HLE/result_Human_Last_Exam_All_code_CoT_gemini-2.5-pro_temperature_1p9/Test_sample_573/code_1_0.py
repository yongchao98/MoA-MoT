import math

def solve_simplices_count():
    """
    Calculates the number of n-simplices for the category N.(Z_N)_k/.

    The number of n-simplices in the nerve of the overcategory Z_{N,k/}
    is the number of non-decreasing sequences of length n+1 from the
    set of objects {k, k+1, ..., N}.
    The size of this set is M = N - k + 1.
    The formula for the number of such sequences is the multiset coefficient
    C(M + (n+1) - 1, n+1), which simplifies to C(N - k + n + 1, n + 1).
    """
    N = 200
    k = 13

    print(f"Given N = {N}, k = {k}:")
    print("The number of n-simplices is given by the formula C(N - k + n + 1, n + 1).")
    print("-" * 30)

    # Store results for the final answer
    results = []

    for n in range(6):
        # Upper argument of the binomial coefficient
        upper_arg = N - k + n + 1
        # Lower argument of the binomial coefficient
        lower_arg = n + 1

        # Calculate the number of simplices
        num_simplices = math.comb(upper_arg, lower_arg)
        results.append(num_simplices)

        # Print the detailed calculation for each n
        print(f"For n = {n}, the number of simplices is:")
        print(f"  C({N} - {k} + {n} + 1, {n} + 1) = C({upper_arg}, {lower_arg}) = {num_simplices}")

if __name__ == '__main__':
    solve_simplices_count()
