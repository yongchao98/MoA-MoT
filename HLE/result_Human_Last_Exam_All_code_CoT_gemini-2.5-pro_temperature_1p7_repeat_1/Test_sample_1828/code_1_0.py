def solve_mad_family_problem():
    """
    This program calculates the difference between the maximal and minimal
    possible cardinality of the set X, based on established theorems in
    ZFC set theory.

    X: Set of cardinalities of uncountable maximal almost disjoint (MAD) families.
    Assumptions: Continuum Hypothesis fails (c > omega_1) and 2^omega_1 = omega_3.
    """

    # Step 1: Determine the minimal possible cardinality of X.
    # In models of Martin's Axiom (MA) with the failure of CH, the set X
    # contains a single element, the continuum cardinality c.
    # It is consistent with the given assumptions to have such a model.
    # For example, a model with MA + c=omega_2 + 2^omega_1=omega_3.
    # In this case, |X| = 1. X cannot be empty.
    min_cardinality_of_X = 1

    # Step 2: Determine the maximal possible cardinality of X.
    # To maximize |X|, we need a model with the largest possible value for c,
    # which is omega_3. A model with c = 2^omega = omega_3 is consistent with
    # the given assumptions.
    # A theorem states it's consistent for MAD families to exist for all
    # cardinals kappa between omega_1 and c with uncountable cofinality.
    # The cardinals in [omega_1, omega_3] are omega_1, omega_2, and omega_3.
    # All three have uncountable cofinality.
    # Therefore, it's consistent that X = {omega_1, omega_2, omega_3}.
    # The size of this set is 3.
    max_cardinality_of_X = 3

    # Step 3: Calculate the difference.
    difference = max_cardinality_of_X - min_cardinality_of_X

    # Print the results, including the final equation.
    print(f"The minimal possible cardinality of X is {min_cardinality_of_X}.")
    print(f"The maximal possible cardinality of X is {max_cardinality_of_X}.")
    print(f"The difference is the result of the equation: {max_cardinality_of_X} - {min_cardinality_of_X} = {difference}.")

# Execute the solution function.
solve_mad_family_problem()