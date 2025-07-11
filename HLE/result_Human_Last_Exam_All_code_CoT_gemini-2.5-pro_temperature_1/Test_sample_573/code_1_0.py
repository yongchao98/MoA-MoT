import math

def solve_simplices_count():
    """
    Calculates the number of n-simplices for the given mathematical construction.
    """
    # Given parameters from the problem description
    N = 200
    k = 13

    print(f"Calculating the number of n-simplices for N={N}, k={k}, and n from 0 to 5.")
    print("-" * 70)

    # Loop for n from 0 to 5
    for n in range(6):
        # An n-simplex in the nerve of the overcategory N.(Z_N)_{k/} corresponds to a
        # non-decreasing sequence of n+1 integers (m_0, ..., m_n) where k <= m_i <= N.
        # The number of such sequences is given by the combination with repetition formula:
        # C(S + r - 1, r), where S is the number of possible integers and r is the length of the sequence.
        # Here, S = N - k + 1 and r = n + 1.
        # The formula is C((N - k + 1) + (n + 1) - 1, n + 1) = C(N - k + n + 1, n + 1).

        # The numbers for the combination formula
        n_total = N - k + n + 1
        k_chosen = n + 1
        
        # Calculate the number of simplices using math.comb
        num_simplices = math.comb(n_total, k_chosen)
        
        # Print the result for each n, showing the numbers in the final equation
        print(f"For n = {n}, the number of {n}-simplices is C({N} - {k} + {n} + 1, {n} + 1) = C({n_total}, {k_chosen}) = {num_simplices}")

solve_simplices_count()