def solve_set_problem():
    """
    Calculates the number of set tuples (S1, S2, S3, S4) satisfying the given conditions.
    """
    # The universe of elements is U = {1, 2, 3, 4, 5}.
    # The conditions are:
    # 1. S1 subset S2 subset S3 subset S4 subset U
    # 2. i in Si for i = 1, 2, 3

    # This problem can be solved by considering the placement of each element
    # of U into one of 5 disjoint regions defined by the set chain:
    # Region 1: S1
    # Region 2: S2 \ S1
    # Region 3: S3 \ S2
    # Region 4: S4 \ S3
    # Region 5: U \ S4

    # Let's determine the number of choices for each element based on the constraints.

    # For element 1: The condition is 1 in S1.
    # This means element 1 must be placed in Region 1.
    choices_for_1 = 1

    # For element 2: The condition is 2 in S2.
    # This means element 2 can be in S1 or in S2 \ S1.
    # So, element 2 can be placed in Region 1 or Region 2.
    choices_for_2 = 2

    # For element 3: The condition is 3 in S3.
    # This means element 3 can be in S1, S2 \ S1, or S3 \ S2.
    # So, element 3 can be placed in Region 1, Region 2, or Region 3.
    choices_for_3 = 3

    # For element 4: There are no specific conditions.
    # It can be placed in any of the 5 regions.
    choices_for_4 = 5

    # For element 5: There are no specific conditions.
    # It can be placed in any of the 5 regions.
    choices_for_5 = 5

    # The total number of ways is the product of the choices for each element.
    total_sets = choices_for_1 * choices_for_2 * choices_for_3 * choices_for_4 * choices_for_5

    print("To find the number of sets (S1, S2, S3, S4), we analyze the possible placements for each element in {1, 2, 3, 4, 5}.")
    print("The choices for each element are determined by the conditions S1 \u2282 S2 \u2282 S3 \u2282 S4 \u2282 {1,2,3,4,5} and i \u2208 Si for i=1,2,3.\n")
    print(f"Number of choices for element 1 (must be in S1): {choices_for_1}")
    print(f"Number of choices for element 2 (must be in S2): {choices_for_2}")
    print(f"Number of choices for element 3 (must be in S3): {choices_for_3}")
    print(f"Number of choices for element 4 (no restrictions): {choices_for_4}")
    print(f"Number of choices for element 5 (no restrictions): {choices_for_5}\n")
    print("The total number of such set tuples is the product of these choices:")
    print(f"Total = {choices_for_1} * {choices_for_2} * {choices_for_3} * {choices_for_4} * {choices_for_5} = {total_sets}")

solve_set_problem()