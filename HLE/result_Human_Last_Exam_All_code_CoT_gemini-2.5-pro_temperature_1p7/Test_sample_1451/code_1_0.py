def solve_trials():
    """
    Calculates the minimum number of trials T(n) for given values of n.

    The formula for T(n) is derived from a recursive strategy and is given by:
    T(n) = 2n - popcount(n)
    where popcount(n) is the number of set bits (1s) in the binary representation of n.
    """

    # The values of n for which T(n) needs to be calculated.
    n_values = [2, 3, 1234, 6712]

    results = []
    for n in n_values:
        # popcount(n) is the number of '1's in the binary representation of n.
        popcount_n = bin(n).count('1')

        # Calculate T(n) using the formula.
        t_n = 2 * n - popcount_n
        results.append(t_n)

    # Print the final list of results, separated by a comma.
    print(",".join(map(str, results)))

solve_trials()