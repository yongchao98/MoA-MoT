def solve_set_theory_problem():
    """
    Calculates the difference between the maximal and minimal possible
    cardinality of the set X based on established results in set theory.
    """

    # The problem asks for the difference between the maximal and minimal possible
    # size of X, where X is the set of cardinalities of uncountable maximal
    # almost disjoint (MAD) families of subsets of omega.
    # The given axioms are ZFC + (2^omega > omega_1) + (2^(omega_1) = omega_3).

    # Minimal possible cardinality of X:
    # Based on consistency results involving Martin's Axiom (MA).
    # It is consistent with the given axioms to have a model of set theory
    # where MA holds. MA forces all MAD families to have the same cardinality, 2^omega.
    # In such a model, the set X contains only one element.
    # Thus, the minimal possible size of X is 1.
    min_card_X = 1

    # Maximal possible cardinality of X:
    # Based on consistency results from Shelah. The assumptions imply that
    # 2^omega must be either omega_2 or omega_3. To maximize the number of
    # possible cardinalities, we consider a model where 2^omega = omega_3.
    # The possible cardinalities for a MAD family are in the interval [omega_1, omega_3].
    # This interval contains three cardinals: omega_1, omega_2, and omega_3.
    # It is consistent with the given axioms to have a model where MAD families
    # of all these three sizes exist.
    # In such a model, the set X has three elements.
    # Thus, the maximal possible size of X is 3.
    max_card_X = 3

    # Calculate the difference.
    difference = max_card_X - min_card_X

    # Output the final equation as requested.
    print(f"{max_card_X} - {min_card_X} = {difference}")

solve_set_theory_problem()