import math

def solve_simplices_count():
    """
    Calculates the number of n-simplices for n=0 to 5
    for the nerve of the over category specified in the problem.
    """
    # Given parameters
    N = 200
    k = 13

    print("Calculating the number of n-simplices for N=200, k=13, and n from 0 to 5.")
    print("-" * 30)

    # Loop through each value of n from 0 to 5
    for n in range(6):
        # The number of n-simplices is given by the combination formula:
        # C(N - k + n + 1, n + 1)
        
        # Upper value for the combination C(m, r)
        m = N - k + n + 1
        # Lower value for the combination C(m, r)
        r = n + 1
        
        # Calculate the number of combinations using math.comb
        result = math.comb(m, r)
        
        # Print the detailed equation and the result for each n
        print(f"For n = {n}:")
        print(f"The number of {n}-simplices is C({N} - {k} + {n} + 1, {n} + 1) = C({m}, {r}) = {result}")
        print("-" * 30)

solve_simplices_count()