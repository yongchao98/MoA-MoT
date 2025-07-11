def solve_cardinality_problem():
    """
    This function calculates the difference between the maximal and minimal
    possible cardinalities of the set X based on set-theoretic principles.
    """

    # The problem specifies that the continuum hypothesis fails and 2^omega_1 = omega_3.
    # Let X be the set of cardinalities of uncountable maximal almost disjoint (MAD) families on omega.

    # Step 1: Determine the value of the continuum, c = 2^omega_0.
    # From set theory, we know (2^omega_0)^omega_1 = 2^(omega_0 * omega_1) = 2^omega_1.
    # So, c^omega_1 = 2^omega_1.
    # By Konig's theorem, c^omega_1 > c.
    # Using the given 2^omega_1 = omega_3, we get omega_3 > c.
    # Since CH fails, c > omega_1, which means c >= omega_2.
    # The only cardinal c satisfying omega_2 <= c < omega_3 is omega_2.
    # Thus, the continuum c = 2^omega_0 = omega_2.
    print("Step 1: Deduced that the continuum c = 2^omega_0 must be omega_2.")

    # Step 2: Determine the possible elements of the set X.
    # The cardinality of any MAD family is at most c. So, any k in X satisfies k <= omega_2.
    # X contains uncountable cardinalities, so k >= omega_1.
    # Therefore, X is a subset of {omega_1, omega_2}.
    print("Step 2: Deduced that X must be a subset of {omega_1, omega_2}.")

    # Step 3: Find the maximal possible cardinality of X.
    # A MAD family of size c = omega_2 always exists. So, omega_2 is always in X.
    # It is consistent with ZFC to have a model where a MAD family of size omega_1 also exists.
    # In such a model, X = {omega_1, omega_2}.
    max_card_X = 2
    print(f"Step 3: The maximal possible cardinality of X is {max_card_X}.")

    # Step 4: Find the minimal possible cardinality of X.
    # As omega_2 is always in X, the minimum size of X is at least 1.
    # It is consistent with ZFC to have a model where there are NO MAD families of size omega_1.
    # In such a model, X = {omega_2}.
    min_card_X = 1
    print(f"Step 4: The minimal possible cardinality of X is {min_card_X}.")

    # Step 5: Calculate the final difference.
    difference = max_card_X - min_card_X
    print("\nFinal Calculation:")
    print(f"The difference between the maximal and minimal possible cardinalities of X is:")
    print(f"{max_card_X} - {min_card_X} = {difference}")

solve_cardinality_problem()