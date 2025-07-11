def solve_limit_problem():
    """
    This function calculates the limit of the expected value of the ratio
    of the number of remaining numbers to n as n approaches infinity.
    The method uses a recurrence relation derived from the problem description.

    Let E_n be the expected number of remaining elements from a row of size n.
    The recurrence relation for E_n is:
        E_n = (2/(n-1)) * sum_{j=0}^{n-2} E_j
    with base cases E_0 = 0 and E_1 = 1.

    We compute the ratio E_n / n for a large n to find the limit.
    """
    from fractions import Fraction

    # We choose a large value for n to get a good approximation of the limit.
    n_limit = 10000

    # E[j] will store the value of E_j
    # S[j] will store the sum S_j = E_0 + ... + E_j
    E = [Fraction(0), Fraction(1)]
    S = [Fraction(0), Fraction(1)]

    # Iteratively compute E_n and S_n up to n_limit
    for n in range(2, n_limit + 1):
        # Calculate E_n using the formula E_n = (2/(n-1)) * S_{n-2}
        E_n = Fraction(2, n - 1) * S[n - 2]
        E.append(E_n)
        
        # Update the sum: S_n = S_{n-1} + E_n
        S_n = S[n - 1] + E_n
        S.append(S_n)

    # The desired ratio is E_n / n for large n.
    ratio = E[n_limit] / n_limit

    print(f"The calculation is based on the recurrence relation E_n = (2/(n-1)) * Sum(E_0 to E_{n-2}).")
    print(f"The ratio E_n / n is computed for n = {n_limit}.")
    print("The limit of the ratio is approximately:")
    print(float(ratio))

solve_limit_problem()