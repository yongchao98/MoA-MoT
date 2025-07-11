import math

def solve_simplices_count():
    """
    Calculates the number of n-simplices for n <= 5, with N=200 and k=13.
    """
    N = 200
    k = 13

    # The objects of the over category are the integers i such that -N <= i <= k.
    # M is the total number of these objects.
    M = k - (-N) + 1

    print(f"For N = {N} and k = {k}, the category of objects is the set of integers from {-N} to {k}.")
    print(f"The total number of objects, M, is {k} - ({-N}) + 1 = {M}.")
    print("The number of n-simplices is the number of non-decreasing sequences of length n+1.")
    print("This is calculated by the formula C(M + n, n + 1).")
    print("-" * 60)

    results = []
    # Loop for n from 0 to 5.
    for n in range(6):
        # The number of n-simplices is given by the binomial coefficient C(M+n, n+1).
        numerator = M + n
        denominator = n + 1
        
        # Calculate the result using math.comb
        result = math.comb(numerator, denominator)
        results.append(result)
        
        # Print the equation and the result for each n.
        print(f"For n = {n}, the number of {n}-simplices is C({M} + {n}, {n} + 1) = C({numerator}, {denominator}) = {result}")

    # The final answer format as requested.
    print(f"\n<<<({', '.join(map(str, results))})>>>")

solve_simplices_count()