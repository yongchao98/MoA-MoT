def solve():
    """
    Calculates T(n), the minimum number of trials needed to decide if we have
    an equal number of real and fake golden bars among 2n bars.
    """
    n_values = [2, 3, 1234, 6712]
    results = []

    for n in n_values:
        # The minimum number of weighings T(n) depends on the parity of n.
        # This is based on the worst-case analysis of a two-stage weighing strategy.
        if n % 2 == 0:
            # For even n, the formula is T(n) = 2n - 1.
            # Example equation for n=2: T(2) = 2 * 2 - 1 = 3
            result = 2 * n - 1
        else:
            # For odd n, the formula is T(n) = 2n - 2.
            # Example equation for n=3: T(3) = 2 * 3 - 2 = 4
            result = 2 * n - 2
        results.append(str(result))

    # Print the final results in the required comma-separated format.
    print(",".join(results))

solve()