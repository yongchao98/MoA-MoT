import math

def solve_topology_problem():
    """
    This function solves the given topology problem by calculating the combination C(k, 2).
    """
    # The problem describes a metric continuum X with k=5 special points.
    # The condition is that no proper subcontinuum of X contains any three of these points.
    # In topology, this identifies X as a "continuum of type A" with respect to these k points.
    # The question asks for the largest number 'n' such that X can be written as an irreducible
    # union of n subcontinua.
    # For a continuum of type A with k points, a known theorem states that this number n
    # is the number of pairs of points, which is the binomial coefficient C(k, 2).

    # Number of special points
    k = 5

    # Number of points in each pairing
    pairs = 2

    # Calculate C(k, 2)
    n = math.comb(k, pairs)

    # To show the full equation as requested, we can also calculate the factorials
    k_factorial = math.factorial(k)
    pairs_factorial = math.factorial(pairs)
    k_minus_pairs_factorial = math.factorial(k - pairs)

    print(f"The number of special points is k = {k}.")
    print(f"The largest number n is the number of pairs of points from this set.")
    print(f"This is calculated by the combination formula C(k, 2).")
    print(f"n = C({k}, {pairs}) = {k}! / ({pairs}! * ({k}-{pairs})!)")
    print(f"n = {k_factorial} / ({pairs_factorial} * {k_minus_pairs_factorial})")
    print(f"n = 120 / (2 * 6)")
    print(f"n = {n}")

solve_topology_problem()