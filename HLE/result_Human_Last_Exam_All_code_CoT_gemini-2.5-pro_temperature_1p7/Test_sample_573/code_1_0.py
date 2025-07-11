import math

def solve_simplices_count():
    """
    Calculates the number of n-simplices for the given mathematical problem.
    """
    N = 200
    k = 13
    
    # An n-simplex in N_.(Z_N)_{k/} corresponds to a sequence of integers
    # (m_0, m_1, ..., m_n) such that k <= m_0 <= m_1 <= ... <= m_n <= N.
    # The set of integers to choose from is {k, k+1, ..., N}.
    # The size of this set is M = N - k + 1.
    M = N - k + 1
    
    print(f"Given N = {N} and k = {k}.")
    print(f"The number of available integers M is N - k + 1 = {N} - {k} + 1 = {M}.")
    print("The number of n-simplices is the number of ways to choose n+1 integers from these M integers with replacement, which is given by the formula C(M+n, n+1).\n")
    
    results = []
    # Loop for n from 0 to 5
    for n in range(6):
        # The number of (n+1)-element multisets from a set of size M is C(M + (n+1) - 1, n+1)
        # which is C(M+n, n+1)
        num_simplices = math.comb(M + n, n + 1)
        results.append(num_simplices)
        
        # Output each number in the final equation as requested.
        print(f"For n = {n}:")
        print(f"  Number of simplices = C({M} + {n}, {n} + 1) = C({M+n}, {n+1}) = {num_simplices}")

    # This part is for the final answer block and not displayed to the user.
    # It prepares the list of numbers for the final answer.
    # print(f"<<<{results}>>>")

solve_simplices_count()