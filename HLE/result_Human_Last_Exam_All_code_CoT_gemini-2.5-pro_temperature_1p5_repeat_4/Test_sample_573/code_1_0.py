import math

def calculate_simplices():
    """
    Calculates the number of n-simplices for n <= 5, with N=200 and k=13,
    for the nerve of the overcategory specified in the problem.

    The problem asks for the number of n-simplices in N.(Z_N)_k/.
    An n-simplex in the nerve of a category corresponds to a sequence of
    n composable morphisms. In our specific overcategory, the objects are integers
    j such that k <= j <= N, and a morphism from j_a to j_b exists if j_a <= j_b.

    Therefore, an n-simplex is a non-decreasing sequence of n+1 integers:
    k <= j_0 <= j_1 <= ... <= j_n <= N.

    The number of such sequences can be found using a stars-and-bars argument.
    It is equivalent to choosing n+1 items from the set of M = N - k + 1 integers
    with replacement. The formula for this is the multiset coefficient:
    C( (n+1) + M - 1, n+1 ) which simplifies to C(n + M, n + 1).
    """
    N = 200
    k = 13

    # The number of objects in the overcategory Z_N/k is M.
    M = N - k + 1

    print(f"For N = {N} and k = {k}, the number of objects M in the overcategory is {N} - {k} + 1 = {M}.")
    print("The number of n-simplices is given by the binomial coefficient C(n + M, n + 1).")
    print("-" * 30)

    results = []
    for n in range(6):
        upper_binom = n + M
        lower_binom = n + 1

        # Use math.comb for accurate calculation of binomial coefficients with large numbers.
        num_simplices = math.comb(upper_binom, lower_binom)
        results.append(num_simplices)

        # Print the equation for each n as requested.
        print(f"For n = {n}, the number of {n}-simplices is C({n} + {M}, {n} + 1) = C({upper_binom}, {lower_binom}) = {num_simplices}")

    return results

if __name__ == '__main__':
    final_results = calculate_simplices()
    # The final answer is the list of results for n=0 to n=5.
    # To conform to the output format, we'll format it as a string.
    # print(f"\n<<<{final_results}>>>")
