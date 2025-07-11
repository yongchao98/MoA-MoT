def solve_cardinality_problem():
    """
    This function calculates the difference between the maximal and minimal
    possible cardinalities of the set X, based on the provided set-theoretic assumptions.

    Let X be the set of cardinalities of uncountable maximal almost disjoint families of subsets of omega.
    Assumptions:
    1. The continuum hypothesis fails (2^aleph_0 > aleph_1).
    2. 2^aleph_1 = aleph_3.
    """

    # As derived from set theory:
    # 1. The minimal possible cardinality of X is 2.
    #    This is because MAD families of size aleph_1 and 2^aleph_0 always exist,
    #    and it is consistent that these are the only possible sizes.
    min_card_X = 2

    # 2. The maximal possible cardinality of X is 3.
    #    This is derived from the fact that 2^aleph_0 must be <= 2^aleph_1 = aleph_3,
    #    so 2^aleph_0 can be at most aleph_3. Under Martin's Axiom, it is consistent
    #    that all cardinals between aleph_1 and 2^aleph_0 are possible sizes for MAD families.
    #    The possible cardinals are aleph_1, aleph_2, and aleph_3.
    max_card_X = 3

    # Calculate the difference
    difference = max_card_X - min_card_X

    # Print the explanation and the final equation
    print(f"Maximal possible cardinality of X: {max_card_X}")
    print(f"Minimal possible cardinality of X: {min_card_X}")
    print(f"The difference is: {max_card_X} - {min_card_X} = {difference}")

solve_cardinality_problem()