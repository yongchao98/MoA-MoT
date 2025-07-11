def solve_weighing_puzzle():
    """
    This function calculates T(n), the minimum number of trials needed to decide
    if we have an equal number of real and fake golden bars among 2n bars.

    The derived formula for T(n) is 2*n - 1. This is based on a worst-case
    analysis of a pairing strategy.

    1.  Perform n weighings by pairing up the 2n bars. This identifies 'u'
        unbalanced pairs (Real, Fake) and 'b' balanced pairs (Real, Real or Fake, Fake).

    2.  The worst-case scenario requires the most additional weighings.
        - If we get at least one unbalanced pair (u > 0), we have a reference
          weight. We then need to check the 'b' balanced pairs, which takes 'b'
          more weighings. Total = n + b = n + (n-u) = 2n - u. The worst case is
          when u is minimal (u=1), giving 2n - 1 weighings.
        - If all pairs are balanced (u = 0), we have n pairs of identical bars.
          We must determine if we have n/2 (Real, Real) and n/2 (Fake, Fake) pairs.
          This subproblem requires an additional n-1 weighings.
          Total = n + (n-1) = 2n - 1.

    Therefore, the minimum number of trials needed in the worst-case is T(n) = 2n - 1.
    """

    def T(n):
        """Calculates the result based on the formula T(n) = 2n - 1."""
        return 2 * n - 1

    values_of_n = [2, 3, 1234, 6712]
    results = []

    print("The formula for the minimum number of trials T(n) is 2*n - 1.")
    print("-" * 35)

    for n in values_of_n:
        result = T(n)
        # As requested, showing each number in the final equation
        print(f"T({n:4}) = 2 * {n:4} - 1 = {result}")
        results.append(str(result))

    print("-" * 35)
    print("The values for T(2), T(3), T(1234), and T(6712) are:")
    print(",".join(results))


solve_weighing_puzzle()

<<<3,5,2467,13423>>>