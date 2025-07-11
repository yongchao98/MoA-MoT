import math

def count_simplices():
    """
    Calculates the number of n-simplices for the given problem.

    An n-simplex in N_.(Z_N)_{k/} corresponds to a sequence of integers
    (c_0, ..., c_n) such that k <= c_0 <= ... <= c_n <= N.

    This is a problem of combinations with repetition. The number of such
    sequences is given by the formula C(m + s - 1, m), where:
    - m = n + 1 is the number of items to choose.
    - s = N - k + 1 is the number of available categories.

    For the given values N=200 and k=13, s = 200 - 13 + 1 = 188.
    The formula becomes C((n + 1) + 188 - 1, n + 1) = C(n + 188, n + 1).
    """
    N = 200
    k = 13
    
    # s is the number of objects in the overcategory Z_N/k
    s = N - k + 1
    
    print(f"For N = {N} and k = {k}, the number of objects in the overcategory is s = {N} - {k} + 1 = {s}.")
    print("The number of n-simplices is calculated using the formula: C(n + s, n + 1)\n")

    # Calculate and print the number of simplices for n from 0 to 5
    for n in range(6):
        # top and bottom values for the binomial coefficient C(top, bottom)
        top = n + s
        bottom = n + 1
        
        result = math.comb(top, bottom)
        
        # Print the full equation for each n
        print(f"For n = {n}:")
        equation = f"Number of {n}-simplices = C({n} + {s}, {n} + 1) = C({top}, {bottom})"
        # The :_ allows for comma as a thousands separator in older python versions
        # The :, is the modern way.
        print(f"{equation} = {result:,}")

count_simplices()