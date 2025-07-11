def solve_set_problem():
    """
    This function calculates the number of set tuples (S1, S2, S3, S4)
    satisfying the given conditions.
    """

    # The problem is to find the number of set tuples (S1, S2, S3, S4) such that:
    # 1. S1 is a subset of S2, S2 of S3, S3 of S4, and S4 of {1, 2, 3, 4, 5}.
    # 2. i is an element of Si for i = 1, 2, 3.

    # We can determine the number of choices for placing each element from the
    # universal set U = {1, 2, 3, 4, 5} into the chain of subsets.
    # The chain creates 5 disjoint regions: S1, S2\S1, S3\S2, S4\S3, U\S4.

    # Number of choices for element 1: must be in S1.
    choices_1 = 1

    # Number of choices for element 2: must be in S2 (i.e., in S1 or S2\S1).
    choices_2 = 2

    # Number of choices for element 3: must be in S3 (i.e., in S1, S2\S1, or S3\S2).
    choices_3 = 3

    # Number of choices for element 4: no constraint, can be in any of the 5 regions.
    choices_4 = 5

    # Number of choices for element 5: no constraint, can be in any of the 5 regions.
    choices_5 = 5

    # The total number of ways is the product of the choices for each element.
    total_sets = choices_1 * choices_2 * choices_3 * choices_4 * choices_5

    # Print the final equation as requested.
    print(f"The total number of sets is the product of the number of choices for each element:")
    print(f"{choices_1} * {choices_2} * {choices_3} * {choices_4} * {choices_5} = {total_sets}")

solve_set_problem()