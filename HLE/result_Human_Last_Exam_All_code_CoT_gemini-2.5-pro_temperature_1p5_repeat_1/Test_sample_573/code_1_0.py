def solve():
    """
    This function calculates the number of n-simplices for the given problem.
    """
    # Parameters from the problem description
    N = 200
    k = 13
    n_max = 5

    # The number of objects in the overcategory is M.
    # These are integers X such that k <= X <= N.
    # The size of this set of objects is M = N - k + 1.
    M = N - k + 1

    # An n-simplex is a non-decreasing sequence (x_0, x_1, ..., x_n)
    # where k <= x_0 <= ... <= x_n <= N.
    # The number of such sequences is given by the multiset coefficient C(M+n, n+1).

    print(f"For N={N}, k={k}, the number of objects is M = {N} - {k} + 1 = {M}.")
    print("The number of n-simplices is calculated using the formula C(M+n, n+1).\n")
    
    # We use an iterative approach to calculate the binomial coefficients
    # to maintain precision with large numbers.
    # S(n) = S(n-1) * (M+n) / (n+1), with S(-1)=1 conceptually, or we can start with S(0).
    
    # Initial value for n=0
    n = 0
    num_simplices = M
    
    # Output the equation for n=0
    # C(M+n, n+1)
    print(f"For n={n}, the number of simplices is C({M}+{n}, {n}+1) = C({M+n}, {n+1}) = {num_simplices}")
    
    # Loop to calculate for n = 1 to n_max
    for n in range(1, n_max + 1):
        # Update the number of simplices using the recurrence relation.
        # The division will be exact.
        num_simplices = num_simplices * (M + n) // (n + 1)
        
        # Output the equation for the current n
        print(f"For n={n}, the number of simplices is C({M}+{n}, {n}+1) = C({M+n}, {n+1}) = {num_simplices}")

solve()