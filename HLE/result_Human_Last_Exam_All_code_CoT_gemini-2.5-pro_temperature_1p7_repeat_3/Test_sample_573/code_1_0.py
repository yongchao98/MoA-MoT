import math

def calculate_simplices():
    """
    Calculates the number of n-simplices for the given problem.

    An n-simplex in N_•(Z_N)_{k/} corresponds to a sequence of integers
    -N <= x_0 <= x_1 <= ... <= x_n <= k.
    This is a "stars and bars" combinatorial problem. The number of such
    sequences is given by the binomial coefficient C(N + k + n + 1, n + 1).
    """
    N = 200
    k = 13

    print(f"Given parameters: N = {N}, k = {k}")
    print("The number of n-simplices is calculated using the formula: C(N + k + n + 1, n + 1)\n")

    for n in range(6):
        # The arguments for the combination C(total, choose)
        total = N + k + n + 1
        choose = n + 1

        # Calculate the result using math.comb for precision with large numbers
        result = math.comb(total, choose)

        # Output the calculation for each value of n
        print(f"For n = {n}:")
        print(f"  Number of {n}-simplices = C({N} + {k} + {n} + 1, {n} + 1)")
        print(f"                       = C({total}, {choose})")
        print(f"                       = {result}\n")

calculate_simplices()