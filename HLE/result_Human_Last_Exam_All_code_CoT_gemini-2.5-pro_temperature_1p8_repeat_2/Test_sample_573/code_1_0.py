import math

def calculate_simplices():
    """
    This function calculates the number of n-simplices for the given problem.

    The number of n-simplices in the specified category corresponds to the number
    of non-decreasing sequences of n+1 integers (x_0, ..., x_n)
    where k <= x_i <= N. This is a classic "stars and bars" problem.
    The formula for the number of such sequences is C(M + n, n + 1),
    where M = N - k + 1 is the number of available integers.
    """
    # Problem parameters
    N = 200
    k = 13

    # The number of available integers for the sequence
    M = N - k + 1

    print(f"For N={N} and k={k}, the number of integers to choose from is M = {N} - {k} + 1 = {M}.")
    print("The number of n-simplices is given by the formula C(M + n, n + 1).\n")
    
    results = []
    # Loop for n from 0 to 5
    for n in range(6):
        
        # Arguments for the combination function C(a, b)
        a = M + n
        b = n + 1
        
        # Calculate the number of n-simplices
        num_simplices = math.comb(a, b)
        results.append(str(num_simplices))
        
        # Print the detailed equation for each n
        print(f"For n={n}, the number of {n}-simplices is "
              f"C({M} + {n}, {n} + 1) = C({a}, {b}) = {num_simplices}")

    # For the final answer format
    final_answer = ",".join(results)
    print(f"\n<<<{final_answer}>>>")

calculate_simplices()
